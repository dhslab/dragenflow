include { DRAGEN_MULTIALIGN    } from '../../modules/local/dragen_multialign.nf'

workflow METHYLATION {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()

    DRAGEN_MULTIALIGN(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_FASTQ_LIST.out.versions)

    emit: 
    ch_versions

}