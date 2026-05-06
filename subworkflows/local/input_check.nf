/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_XLSX_TO_CSV } from '../../modules/local/convert_xlsx_to_csv'

/*
========================================================================================
    SUBWORKFLOW TO CHECK INPUTS
========================================================================================
*/

workflow INPUT_CHECK {

    take:
    input          //  string: Path to input samplesheet
    ch_fastq_list  //  channel: [ path(file) ]

    main:
    ch_samples  = Channel.empty()
    ch_versions = Channel.empty()

    /*
    ================================================================================
                        Process input MGI samplesheet
    ================================================================================
    */

    if (input) {
        // Verify input samplesheet has a file extension in [xlsx,csv,tsv]
        if (hasExtension(input, 'xlsx')) {
            CONVERT_XLSX_TO_CSV (
                Channel.fromPath(input, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(CONVERT_XLSX_TO_CSV.out.versions)

            ch_input = CONVERT_XLSX_TO_CSV.out.csv

        } else if (hasExtension(input, 'csv') || hasExtension(input, 'tsv')) {
            ch_input = Channel.fromPath(input, checkIfExists: true)
        } else {
            error("Input samplesheet input does not end in `.{xlsx,csv,tsv}`!")
        }
    } else {
        ch_input = Channel.empty()
    }

    /*
    ================================================================================
                        Process input FastQ list
    ================================================================================
    */

    ch_samples = ch_samples.mix(
                    ch_fastq_list
                        .filter{ it != [] }
                        .flatMap{
                            def data = parseFastqList(it)
                            def requiredColumns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
                            data.collect{
                                if (!it.keySet().containsAll(requiredColumns)) {
                                    error("Missing required columns in input FastQ list!")
                                }

                                def R1 = file(it['Read1File'], checkIfExists: true)
                                def R2 = file(it['Read2File'], checkIfExists: true)

                                def regexPattern = /\w\d{2}-\d+/
                                def meta = [
                                    id  : it.RGSM,
                                    acc : it.RGSM.findAll(regexPattern)[0] ?: it.RGSM,
                                    RGSM: it.RGSM
                                ]

                                [ meta.acc, meta, [ R1, R2 ] ]
                            }
                        }
                        .groupTuple()
                        .combine(
                            ch_fastq_list
                                .filter{ it != [] }
                                .map{
                                    def data = parseFastqList(it)
                                    data.each{
                                        if (it) {
                                            it['Read1File'] = "fastq_files/${it['Read1File'].split('/')[-1]}"
                                            it['Read2File'] = "fastq_files/${it['Read2File'].split('/')[-1]}"
                                        }
                                    }

                                    if (data) {
                                        def header = data[0].keySet().join(',')
                                        def content = data.collect { it.values().join(',') }.join('\n')

                                        [ [ header ], [ content ] ]
                                    } else {
                                        [ [], [] ]
                                    }
                                }
                                .flatten()
                                .collectFile(
                                    name   : "updated_fastq_list.csv",
                                    newLine: true,
                                    sort   : 'index'
                                )
                        )
                        .map{ id, meta, reads, fastq_list -> [ meta.first(), reads.flatten(), fastq_list, [] ] }
                )

    emit:
    samples  = ch_samples   // channel: [ val(sample_info), path(reads), path(fastq_list), path(empty list - placeholder for alignment file) ]
    input    = ch_input     // channel: [ path(file) ]
    versions = ch_versions  // channel: [ path(file) ]

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Get file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
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
