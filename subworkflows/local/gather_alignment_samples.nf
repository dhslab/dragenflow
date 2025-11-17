/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARTITION_ALIGNMENT_FILE         } from '../../modules/local/partition_alignment_file'
include { CONVERT_ALIGNMENT_FILE_TO_FASTQ  } from '../../modules/local/convert_alignment_file_to_fastq'
include { CREATE_FASTQ_LIST                } from '../../modules/local/create_fastq_list'


// Parse CSV files from 'ch_samples_to_align' for each sample into 'meta'
def generateMetaFromCsv(csv_file) {
    def lines = csv_file.text.readLines()
    def headers = lines[0].split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')

    return lines.drop(1).collect { line ->
        def fields = line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
        [headers, fields].transpose().collectEntries { k, v -> v ? [(k): v] : [:] }
    }.findAll { it }
}

// Parse FastQ list
def parseFastqList(file) {
    def separator = file.toString().endsWith("tsv") ? '\t' : ','
    def lines = file.readLines()
    def headers = lines.first().split(separator)
    lines.drop(1).collect{ line ->
        [headers, line.split(separator)].transpose().collectEntries{ it }
    }
}

/*
========================================================================================
    SUBWORKFLOW TO GATHER ALIGNMENT SAMPLES
========================================================================================
*/

workflow GATHER_ALIGNMENT_SAMPLES {

    take:
    ch_samples_to_align  // channel: [ path(file) ]
    ch_demux_path        // channel: [ val(meta), path(demux_path) ]
    ch_cram_reference    // channel: [ path(file) ]

    main:
    ch_versions          = Channel.empty()
    ch_gathered_fastqs   = Channel.empty() // channel: [ val(meta), path(read1), path(read2), path(runinfo.xml) ]

    // This is the output of this subworkflow and is the format for DRAGEN alignment
    ch_gathered_samples  = Channel.empty() // [ val(sample_info), path(reads), path(fastq_list), path(empty list - placeholder for alignment file) ]

    // Channel of meta data for alignment samples
    ch_sample_alignment_meta = ch_samples_to_align
                                .map{ generateMetaFromCsv(it) }
                                .flatten()

    //
    // Collect reads, fastq_list, and runinfo.
    //
    ch_gathered_fastqs = ch_gathered_fastqs
        .mix (
            ch_sample_alignment_meta
                .filter{ it.read1 && it.read2 }
                .map{ meta -> 
                    def newMeta = meta.subMap('id', 'individual_id','sample_type','sample_id','sex')
                    [ newMeta, file(meta.read1, checkIfExists: true), file(meta.read2, checkIfExists: true), [] ] 
                },
            ch_sample_alignment_meta
                .filter{ it.fastq_list }
                .flatMap{ meta -> 
                    def newMeta = meta.subMap('id', 'individual_id','sample_type','sample_id','sex')
                    def requiredColumns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
                    def fastq_list = file(meta.fastq_list, checkIfExists: true)                    
                    def data = parseFastqList(fastq_list)
                    data = data.findAll{ it.RGSM == newMeta.id }                    
                    data.collect{
                        if (!it.keySet().containsAll(requiredColumns)) {
                            error("Missing required columns in input FastQ list!")
                        }
                        def R1 = file(it['Read1File'], checkIfExists: true)
                        def R2 = file(it['Read2File'], checkIfExists: true)                    
                        [ newMeta, R1, R2, [] ]
                    }
                }
                .filter{ it!= [] },
            ch_demux_path
                .mix(
                    ch_sample_alignment_meta
                            .filter{ it.demux_path }
                            .map { meta -> 
                                def newMeta = meta.subMap('id', 'individual_id','sample_type','sample_id','sex')
                                [ newMeta, file(meta.demux_path, checkIfExists: true) ] 
                            }
                )
                .filter{ it != [] }
                .flatMap{ meta, demux_path -> 
                    def requiredColumns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
                    def runinfo = file("${demux_path}/Reports/RunParameters.xml").with { it.exists() ? it : [] }
                    def fastq_list = file("${demux_path}/Reports/fastq_list.csv", checkIfExists: true)
                    def data = parseFastqList(fastq_list)
                    // extract data items where RGSM == meta.id
                    data = data.findAll{ it.RGSM == meta.id }
                    data.collect{
                        if (!it.keySet().containsAll(requiredColumns)) {
                            error("Missing required columns in input FastQ list!")
                        }
                        def R1 = file("${demux_path}" + "/" + it['Read1File'].split('/')[-1], checkIfExists: true)
                        def R2 = file("${demux_path}" + "/" + it['Read2File'].split('/')[-1], checkIfExists: true)
                        [ meta, R1, R2, runinfo ]
                    }
                }
                .filter{ it!= [] }
        )
    
    //
    // MODULE: Create fastq_list with local/staged fastq paths and appropriate metadata
    //
    CREATE_FASTQ_LIST (
        ch_gathered_fastqs
    )
    ch_versions = ch_versions.mix(CREATE_FASTQ_LIST.out.versions)

    //
    // SUBWORKFLOW: Verify and parse fastq_list files and staged reads
    //
    ch_gathered_samples = ch_gathered_samples
        .mix(
            ch_gathered_fastqs
                .filter{ it != [] }
                .map{ meta, read1, read2, runinfo -> [ meta, [ read1, read2 ] ] }
                .groupTuple()
                .map { meta, reads -> [ meta, reads.flatten() ] }
            .join(
                CREATE_FASTQ_LIST.out.fastq_list
                .map { meta, fastq_list -> [ meta.id, meta ] }
                .join(
                    CREATE_FASTQ_LIST.out.fastq_list  // NOTE (DHS): This creates a fastq_list will all samples in this run that is reused.
                        .map{ meta, fastq_list ->
                            def data = parseFastqList(fastq_list)
                            data.each{
                                if (it) {
                                    it['Read1File'] = "fastq_files/${it['Read1File'].split('/')[-1]}"
                                    it['Read2File'] = "fastq_files/${it['Read2File'].split('/')[-1]}"
                                }
                            }
                            if (data) {
                                def header = data[0].keySet().join(',')
                                def content = data.collect { it.values().join(',') }.join('\n')

                                [ meta, header + '\n' + content ]
                            } else {
                                []
                            }
                        }
                        .collectFile(
                            { meta, content -> [ "${meta.id}.fastq_list.csv", content + '\n' ] },
                            keepHeader: true,
                            newLine: false
                        )
                        .map{ filename -> [ filename.getName().replaceAll("\\.fastq_list\\.csv", ""), filename ] }
                )
                .map { id, meta, fastq_list -> [ meta, fastq_list ] }
            )
            .map{ meta, reads, fastq_list -> [ meta, reads.flatten(), fastq_list, [] ] }
        )

    /*
    //
    // Get CRAM/BAM files that can be directly realigned
    //
    ch_ = ch_gathered_samples.mix(
                            ch_sample_alignment_meta
                                .map { meta -> 
                                    def validKeys = ['read1', 'read2', 'demux_path', 'cram', 'bam', 'fastq_list']
                                    def item_count = meta.count { key, value -> (key in validKeys) && value }
                                    if (item_count == 1){
                                        return meta
                                    }
                                }
                                .filter{ (it?.bam || it?.cram) && !(it?.read1 || it?.read2 || it?.demux_path || it?.fastq_list) }
                                .map{ [ it, [], [], file(it.bam ? "${it.bam}*" : "${it.cram}*") ] }
                            )
   
    //
    // SUBWORKFLOW: Partition BAM/CRAM files that need to be converted to fastq and realigned. 
    //              This retains run information and can handle bam/cram for the same sample.
    //
    ch_crams_to_partition = ch_sample_alignment_meta
            .map { meta -> 
                def validKeys = ['read1', 'read2', 'demux_path', 'cram', 'bam', 'fastq_list']
                def item_count = meta.count { key, value -> (key in validKeys) && value }
                if (item_count > 1){
                    return meta
                }
            }
            .filter{ it.bam || it.cram }
            .flatMap { meta ->
                    def outputs = []
                    if (meta.cram) {
                        outputs.add([meta, file("${meta.cram}*")])
                    }
                    if (meta.bam) {
                        outputs.add([meta, file("${meta.bam}*")])
                    }
                    return outputs
            }

    PARTITION_ALIGNMENT_FILE (
        ch_crams_to_partition,
        ch_cram_reference
    )
    ch_versions = ch_versions.mix(PARTITION_ALIGNMENT_FILE.out.versions)
    
    //
    // MODULE: Convert CRAM files to FastQ
    //
    CONVERT_ALIGNMENT_FILE_TO_FASTQ (
        PARTITION_ALIGNMENT_FILE.out.cram_files
            .flatMap { meta, values ->
                values.collect { val -> [meta, val] }
            },
        ch_cram_reference
    )
    ch_versions = ch_versions.mix(CONVERT_ALIGNMENT_FILE_TO_FASTQ.out.versions)

    ch_gathered_fastqs = ch_gathered_fastqs.mix (
        CONVERT_ALIGNMENT_FILE_TO_FASTQ.out.fastqs
            .map{ meta, read1, read2 -> [ meta, read1, read2, [] ] }
    )
    */

    //ch_gathered_samples.dump(tag: 'gather_alignment_samples:ch_gathered_samples', pretty: true)

    emit:
    samples  = ch_gathered_samples  // channel: [ val(meta), path(reads), path(fastq_list), path(alignment_file) ]
    versions = ch_versions           // channel: [ path(file) ]

}
