include { DRAGEN_MULTIALIGN as DRAGEN_SOMATIC } from '../../modules/local/dragen_multialign.nf'

workflow SOMATIC {
    take:
    input_data
    dragen_inputs
    
    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_SOMATIC(input_data, dragen_inputs)
    ch_dragen_output = ch_dragen_output.mix(DRAGEN_SOMATIC.out.dragen_output)
    ch_versions = ch_versions.mix(DRAGEN_SOMATIC.out.versions)

    emit:
    dragen_output = ch_dragen_output
    versions = ch_versions

}
