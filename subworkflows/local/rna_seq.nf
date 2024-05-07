include { DRAGEN_MULTIALIGN as DRAGEN_RNA } from '../../modules/local/dragen_multialign.nf'
include { ANNOTATE_RNASEQ                 } from '../../modules/local/annotate_rnaseq.nf'

workflow RNASEQ {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_RNA(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_RNA.out.versions)

    ANNOTATE_RNASEQ(DRAGEN_RNA.out.dragen_output, dragen_inputs)
    ch_versions = ch_versions.mix(ANNOTATE_RNASEQ.out.versions)
    
    emit:
    dragen_output = DRAGEN_RNA.out.dragen_output
    versions = ch_versions
}