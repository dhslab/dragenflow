/*

This subworkflow accepts a channel of sample metadata and inputs and returns a channel of FASTQ files and their associated metadata.
The sample metadata and files must be in the following format:

If cram/bam inputs are given when flowcell/lane/indexes are also given then the CRAM/BAM files converted to FASTQ files
so these can be combined with new data from the DEMUX subworkflow.

The output channel contains the following format:

The rgpl has the following format: 
    <instrument>.<side>.<flowcellid>.<lane>.<flowcelltype>.<reagentlot>.<runrecipe>.<i7index>.<i5index>

Note the some fields can be empty and will be if the input metadata is missing the relevant information.

*/

include { MAKE_FASTQLIST as MAKE_FASTQLIST_FROM_PATH } from '../../modules/local/make_fastqlist.nf'
include { MAKE_FASTQLIST as MAKE_FASTQLIST_FROM_LIST } from '../../modules/local/make_fastqlist.nf'
include { MAKE_FASTQLIST as MAKE_FASTQLIST_FROM_READS } from '../../modules/local/make_fastqlist.nf'
include { MAKE_FASTQLIST as MAKE_FASTQLIST_FROM_CRAMS } from '../../modules/local/make_fastqlist.nf'

def parse_fastqlist(LinkedHashMap row) {
       
    def read1 = new File(row.Read1File)
    def read2 = new File(row.Read2File)

    if (!file(read1).exists()) { 
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${read1}"
    }
    if (!file(read2).exists()) { 
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${read2}"
    }

    row.Read1File = read1
    row.Read2File = read2

    fastq_meta = [ row.RGID, row.RGSM, row.RGLB, row.Lane, row.RGPU, file(read1), file(read2) ]

    return row

}

def parseCSV(filePath) {
    List<Map<String, String>> data = []

    // Read the file and split each line
    filePath.withReader { reader ->
        // Read the header line to use as keys for the map
        def headers = reader.readLine().split(',')

        // Process each subsequent line
        reader.splitEachLine(',') { values ->
            def row = [:]
            headers.eachWithIndex { header, i ->
                row[header.trim()] = values[i].trim()
            }
            data.add(row)
        }
    }

    return data
}

workflow GATHER_FASTQS {
    take:
    ch_sample_meta

    main:

    ch_fastq_lists = Channel.empty()
    ch_versions = Channel.empty()

    //
    // Get FASTQs from demux_path
    // 
    ch_sample_meta
    .map { meta, inputs -> 
        if (inputs.demux_path) {
            [ meta, inputs.demux_path ]
        }
    }
    .transpose()
    .set { ch_demux_path }

    MAKE_FASTQLIST_FROM_PATH(ch_demux_path, "demux_path")
    ch_versions = ch_versions.mix(MAKE_FASTQLIST_FROM_PATH.out.versions)

    MAKE_FASTQLIST_FROM_PATH.out.fastqlist
    .collectFile(keepHeader:false,newLine:false,skip:1) { meta, fastqlist -> 
        [ "${meta.id}.from_path_fastq_list.csv", fastqlist ]
    }
    .map { fqfile -> [ fqfile.getName().toString().split('\\.')[0], fqfile ] }
    .join(ch_sample_meta.map { meta, inputs -> [ meta.id, meta ] })
    .map { id, fastqlist, meta -> 
        [ meta, fastqlist ]
    }
    .mix(ch_fastq_lists)
    .set { ch_fastq_lists }

    // Get FASTQs from fastq_list
    //
    ch_sample_meta
    .map { meta, inputs -> 
        if (inputs.fastq_list) {
            [ meta, inputs.fastq_list ]
        }
    }
    .transpose()
    .set { ch_fastq_list }

    MAKE_FASTQLIST_FROM_LIST(ch_fastq_list, "fastq_list")
    ch_versions = ch_versions.mix(MAKE_FASTQLIST_FROM_LIST.out.versions)
    
    MAKE_FASTQLIST_FROM_LIST.out.fastqlist
    .collectFile(keepHeader:false,newLine:false,skip:1) { meta, fastqlist -> 
        [ "${meta.id}.from_list_fastq_list.csv", fastqlist ]
    }
    .map { fqfile -> [ fqfile.getName().toString().split('\\.')[0], fqfile ] }
    .join(ch_sample_meta.map { meta, inputs -> [ meta.id, meta ] })
    .map { id, fastqlist, meta -> 
        [ meta, fastqlist ]
    }
    .mix(ch_fastq_lists)
    .set { ch_fastq_lists }


    //
    // Get FASTQLIST from reads
    //
    ch_sample_meta
    .map { meta, inputs -> 
        if (inputs.reads) {
            [ meta, inputs.reads ]
        }
    }
    .transpose()
    .map { meta, reads -> 
        [ meta, [reads.read1, reads.read2] ]
    }
    .set { ch_reads }

    MAKE_FASTQLIST_FROM_READS(ch_reads, "reads")
    ch_versions = ch_versions.mix(MAKE_FASTQLIST_FROM_READS.out.versions)

    MAKE_FASTQLIST_FROM_READS.out.fastqlist
    .collectFile(keepHeader:false,newLine:false,skip:1) { meta, fastqlist -> 
        [ "${meta.id}.from_reads_fastq_list.csv", fastqlist ]
    }
    .map { fqfile -> [ fqfile.getName().toString().split('\\.')[0], fqfile ] }
    .join(ch_sample_meta.map { meta, inputs -> [ meta.id, meta ] })
    .map { id, fastqlist, meta -> 
        [ meta, fastqlist ]
    }
    .mix(ch_fastq_lists)
    .set { ch_fastq_lists }

    //
    // Get CRAMS and BAMS if there are reads/fastqs for the same sample 
    //
    ch_sample_meta
    .map { meta, inputs -> 
        if (inputs.cram){
            [ meta, inputs.cram ]
        }
    }
    .transpose()
    .set { ch_crams }

    //
    // SPLIT CRAMS, CONVERT TO FASTQ AND CREATE FASTQLIST CHANNEL
    //
    SPLIT_CRAM(ch_crams, Channel.value(params.fasta))
    ch_versions = ch_versions.mix(SPLIT_CRAM.out.versions)
    
    FASTQ_FROM_CRAM(SPLIT_CRAM.out.crams.transpose(), Channel.value(params.fasta))
    ch_versions = ch_versions.mix(FASTQ_FROM_CRAM.out.versions)

    MAKE_FASTQLIST_FROM_CRAMS(FASTQ_FROM_CRAM.out.fastqs, "reads")
    ch_versions = ch_versions.mix(MAKE_FASTQLIST_FROM_READS.out.versions)

    MAKE_FASTQLIST_FROM_CRAMS.out.fastqlist
    .collectFile(keepHeader:false,newLine:false,skip:1) { meta, fastqlist -> 
        [ "${meta.id}.from_crams_fastq_list.csv", fastqlist ]
    }
    .map { fqfile -> [ fqfile.getName().toString().split('\\.')[0], fqfile ] }
    .join(ch_sample_meta.map { meta, inputs -> [ meta.id, meta ] })
    .map { id, fastqlist, meta -> 
        [ meta, fastqlist ]
    }
    .mix(ch_fastq_lists)
    .set { ch_fastq_lists }

emit:
    fastqs = ch_fastq_lists
    versions = ch_versions
}


//
// Processes
//

process SPLIT_CRAM {
    tag "$meta.id"
    label 'process_medium_multithread'
    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
        tuple val(meta), path(cramfile)
        path(reference)

    output:
        tuple val(meta), path("*.cram"),    emit: crams
        path("versions.yml"),               emit: versions

    script:
        """
        samtools split -@ ${task.cpus} -f '%!.cram' --write-index --reference ${reference} ${cramfile}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        """

    stub:
        """
        cp ${projectDir}/assets/stub/split_crams/${meta.id}/*.cram .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        END_VERSIONS
        """
}

process FASTQ_FROM_CRAM {
    tag "$meta.id"
    label 'process_high_multithread'
    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
        tuple val(meta), path(cramfile)
        path(reference)

    output:
        tuple val(meta), path("*{R1,R2}.fastq.gz", arity: 2), emit: fastqs
        path("versions.yml"),    emit: versions

    script:
    def threads = (task.cpus / 2) - 1

    """
    samtools collate -@ ${threads} --reference ${reference} -u -O ${cramfile} | \
    samtools fastq -@ ${threads} -1 \$(basename ${cramfile} .cram).R1.fastq.gz -2 \$(basename ${cramfile} .cram).R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    cp ${projectDir}/assets/stub/cram_to_fastq/${meta.id}/*.fastq.gz .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

}
