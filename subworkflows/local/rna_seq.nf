include { DRAGEN_MULTIALIGN } from '../../modules/local/dragen_multialign.nf'
include { ANNOTATE_RNASEQ   } from '../../modules/local/annotate_rnaseq.nf'

workflow RNASEQ {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_MULTIALIGN(input_data, dragen_inputs)
    ch_dragen_output = ch_dragen_output.mix(DRAGEN_MULTIALIGN.out.dragen_output)
    ch_versions = ch_versions.mix(DRAGEN_MULTIALIGN.out.versions)

    ANNOTATE_RNASEQ(ch_dragen_output, dragen_inputs)

    emit:
    ch_versions
}