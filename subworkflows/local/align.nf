include { DRAGEN_MULTIALIGN as DRAGEN_ALIGN     } from '../../modules/local/dragen_multialign.nf'

workflow ALIGN {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_ALIGN(input_data, dragen_inputs)
    ch_dragen_output = ch_dragen_output.mix(DRAGEN_ALIGN.out.dragen_output)
    ch_versions = ch_versions.mix(DRAGEN_ALIGN.out.versions)

    emit:
    dragen_output = ch_dragen_output
    versions = ch_versions
}
