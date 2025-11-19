
workflow PREPARE_SOMATIC_FASTQS {
    take:
    ch_prepare_somatic_fastqs_input  //  channel: [ id, [ path(fastqlist1), fastqlist2, ... ] ]

    main:
    ch_prepare_somatic_fastqs_output = Channel.empty()

    ch_prepare_somatic_fastqs_input
        .branch {
            tumor: it[0].sample_type == 'tumor'
            normal: it[0].sample_type == 'normal'
        }.set { ch_prepare_somatic_fastqs_samples }

    ch_prepare_somatic_fastqs = ch_prepare_somatic_fastqs_samples.tumor // first get tumor samples
        .map { meta, reads, fastqlist, alignment_files ->
            [ meta.individual_id, meta, reads, fastqlist ]
        }
        .cross( // and join with normal samples
            ch_prepare_somatic_fastqs_samples.normal
            .map { meta, reads, fastqlist, alignment_files ->
                [ meta.individual_id, meta, reads, fastqlist ]
            }
        ).map { tumor, normal -> [ tumor[1], tumor[2], tumor[3], normal[1], normal[2], normal[3] ] }
        // on joined samples, set tumor and normal id to meta and combine reads.
        .map { tumor_meta, tumor_reads, tumor_fastqlist, normal_meta, normal_reads, normal_fastqlist ->
            def new_meta = [:]
            new_meta['id'] = tumor_meta['id']
            new_meta['sex'] = tumor_meta['sex']
            new_meta['tumor_id'] = tumor_meta['id']
            new_meta['normal_id'] = normal_meta['id']
            reads = tumor_reads + normal_reads
            [ new_meta, reads, [ tumor_fastqlist, normal_fastqlist ] ]
        }
        // now separate reads and fastqlists
        .multiMap { meta, reads, fastqlists -> 
            reads: [ meta, reads ]
            fastqlist: [ meta.id, meta, fastqlists ]
        }
        
    ch_prepare_somatic_fastqs_output = ch_prepare_somatic_fastqs.reads
        .join(CAT_FASTQLISTS(ch_prepare_somatic_fastqs.fastqlist).fastqlist)
        .map { meta, reads, fastqlist ->
            [ meta, reads, fastqlist, [] ]
        }

    emit: 
    samples = ch_prepare_somatic_fastqs_output

}

workflow CAT_FASTQLISTS {
    take:
    ch_cat_fastqlists_input  //  channel: [ id, [ path(fastqlist1), fastqlist2, ... ] ]

    main:
    ch_cat_fastqlists_output = Channel.empty()

    ch_cat_fastqlists_input
        .map { id, meta, fastqlists -> [ id, meta ] }
        .join(
            ch_cat_fastqlists_input
            .map { id, meta, fastqlists ->
                def data = fastqlists.collectMany { file -> parseFastqList(file) }
                if (!data.isEmpty()) {
                    def header = data[0].keySet().join(',')
                    def content = data.collect { row ->
                        row.values().join(',')
                    }.join('\n')
                    [ id, header + '\n' + content + '\n' ]
                }
            }
            .collectFile(
                { id, content -> [ "${id}.fastq_list.csv", content + '\n' ] },
                keepHeader: true,
                newLine: false
            )
            .map{ filename -> [ filename.getName().replaceAll("\\.fastq_list\\.csv", ""), filename ] }
        )
        .map { id, meta, fastqfile -> [ meta, fastqfile ] }
        .set { ch_cat_fastqlists_output }

    emit:
    fastqlist = ch_cat_fastqlists_output

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Parse FastQ list
def parseFastqList(file) {
    def separator = file.toString().endsWith("tsv") ? '\t' : ','
    def lines = file.readLines()
    def headers = lines.first().split(separator)
    lines.drop(1).collect{ line ->
        [headers, line.split(separator)].transpose().collectEntries{ it }
    }
}