include { SAMPLESHEET_CHECK              } from '../../modules/local/samplesheet_check.nf'
include { MAKE_FASTQLIST                 } from '../../modules/local/make_fastqlist.nf'

workflow INPUT_CHECK {
    take:
    master_samplesheet
    data_path

    main:

    ch_mastersheet        = Channel.empty()
    ch_input_data         = Channel.empty()
    ch_dragen_outputs     = Channel.empty()
    ch_versions           = Channel.empty()

    // Runs a python script that parses the sample sheet and adds key metadata, 
    // including index sequences, flowcell, and lane. If fastq_list.csv files are passed,
    // these are parsed also and read1/read2 pairs are returned. The output is then channelified.
    SAMPLESHEET_CHECK ( master_samplesheet, data_path )
    .csv
    .splitCsv ( header:true, sep:',' )
    .map { create_master_samplesheet(it) }
    .set { ch_mastersheet }

    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    // Organize reads into a fastq list string (to be written to a file) and read1/read2 pairs.
    ch_mastersheet
    .map { meta -> 
        if (meta.read1 != null && meta.read2 != null){
            def new_meta = meta.subMap('id', 'uid','sex','sample_type','sample_id','assay')
            def rgid = meta.flowcell + '.' + meta.i7index + '.' + meta.i5index + '.' + meta.lane 
            def rglb = meta.id + '.' + meta.i7index + '.' + meta.i5index
            [ new_meta, [ rgid, meta.id, rglb, meta.lane, file(meta.read1), file(meta.read2) ] ]
        }
    }
    .groupTuple()
    .map { meta, fqlist -> 
            def fileList = ['RGID,RGSM,RGLB,Lane,Read1File,Read2File']
            def read1 = []
            def read2 = []

            // Create data rows
            for (int i = 0; i < fqlist.size(); i++) {
                def row = fqlist[i]
                read1 << file(row[4])
                read2 << file(row[5])
                fileList << [ row[0], row[1], row[2], row[3], row[4].toString().split('/')[-1], row[5].toString().split('/')[-1] ].join(',')
            }
            return [ meta, fileList.join('\n'), read1, read2 ]
    }
    .set { ch_fastqs }

    ch_fastqs
    .collectFile{ meta, fqlist, read1, read2 -> 
        [ "${meta.id}.fastq_list.csv", fqlist + '\n' ]
    }
    .map { fqfile -> 
        [ fqfile.getName().toString().split('\\.')[0], file(fqfile) ]
    }
    .set { ch_fastq_list_files }

    ch_fastqs
    .map { meta, fqlist, read1, read2 ->
        [ meta.id, meta ]
    }
    .join(ch_fastq_list_files)
    .map { id, meta, fqlist -> 
        [ meta, fqlist ]
    }
    .set { ch_fastq_lists }

    // Put read1 and read2 files into separate channels.
    ch_fastqs
    .map { meta, fqlist, read1, read2 -> [ meta, read1 ] }
    .transpose()
    .set { ch_read1 }

    ch_fastqs
    .map { meta, fqlist, read1, read2 -> [ meta, read2 ] }
    .transpose()
    .set { ch_read2 }

    ch_fastq_lists
    .concat(ch_read1,ch_read2)
    .groupTuple()
    .map { meta, files -> 
        [ meta, "fastq", files ] 
    }
    .set { ch_fastq_list }

    ch_input_data = ch_input_data.mix(ch_fastq_list)

    // Organize cram files into a channel.
    ch_mastersheet
    .map { meta -> 
        if (meta.cram != null){
            def new_meta = meta.subMap('id', 'uid','sex','sample_type','sample_id','assay')
            new_meta.cram = file(meta.cram).getName()
            [ new_meta, 'cram', [ file(meta.cram), file(meta.cram + '.crai')  ] ]
        }
    }
    .set { ch_cram }

    ch_input_data = ch_input_data.mix(ch_cram)

    // Organize bam files into a channel.
    ch_mastersheet
    .map { meta -> 
        if (meta.bam != null){
            def new_meta = meta.subMap('id', 'uid','sex','sample_type','sample_id','assay')
            new_meta.bam = file(meta.bam).getName()
            [ new_meta, 'bam', [ file(meta.bam), file(meta.bam + '.bai') ] ]
        }
    }
    .set { ch_bam }

    ch_input_data = ch_input_data.mix(ch_bam)

    emit:
    input_data = ch_input_data
    versions = ch_versions

}

def create_master_samplesheet(LinkedHashMap row) {

    def meta = [:]
    meta.id             = row.id
    meta.uid            = row.uid ?: null
    meta.sex            = row.sex ?: null
    meta.sample_type    = row.sample_type ?: null
    meta.sample_id      = row.sample_id ?: null
    meta.assay          = row.assay ?: null
    meta.i7index        = row.i7index ?: null
    meta.i5index        = row.i5index ?: null
    meta.flowcell       = row.flowcell ?: null
    meta.lane           = row.lane ?: null

    if (row.read1) {
        if (!file(row.read1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 file does not exist!\n${row.read1}"
        }
        meta.read1 = file(row.read1)
    }
    if (row.read2) {
        if (!file(row.read2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 file does not exist!\n${row.read2}"
        }
        meta.read2 = file(row.read2)
    }

    if (row.cram) {
        if (!file(row.cram).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Cram file does not exist!\n${row.cram}"
        }
        meta.cram = file(row.cram)
    }

    if (row.bam) {
        if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
        }
        meta.bam = file(row.bam)
    }
    return meta
}
