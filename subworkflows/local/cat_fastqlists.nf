
workflow CAT_FASTQLISTS {
    take:
    ch_fastqlists2cat  //  channel: [ id, [ path(fastqlist1), fastqlist2, ... ] ]

    main:
    ch_catfastqlist_output = Channel.empty()
    ch_fastqlists2cat
        .map { id, meta, fastqlists -> [ id, meta ] }
        .join(
            ch_fastqlists2cat
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
        .set { ch_catfastqlist_output }

    emit:
    fastqlist = ch_catfastqlist_output

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