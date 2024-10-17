include { DRAGEN_MULTIALIGN as DRAGEN_METHYLATION    } from '../../modules/local/dragen_multialign.nf'

workflow METHYLATION {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()

    DRAGEN_METHYLATION(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_METHYLATION.out.versions)

    // Need to add a process to convert dragen methylation output to a bed file.
    
    emit: 
    dragen_outputs = DRAGEN_METHYLATION.out.dragen_output
    versions = ch_versions

}