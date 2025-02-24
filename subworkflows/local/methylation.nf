include { DRAGEN_MULTIALIGN as DRAGEN_METHYLATION    } from '../../modules/local/dragen_multialign.nf'
include { MAKE_METH_BED } from '../../modules/local/make_meth_bed.nf'
include { MAKE_METH_BIGWIG } from '../../modules/local/make_meth_bigwig.nf'

workflow METHYLATION {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()

    DRAGEN_METHYLATION(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_METHYLATION.out.versions)

    // Need to add a process to convert dragen methylation output to a bed file.
    MAKE_METH_BED(DRAGEN_METHYLATION.out.dragen_output)
    ch_versions = ch_versions.mix(MAKE_METH_BED.out.versions)

    // make channel for reference and index
    ch_fasta_reference = params.fasta
    ? Channel.fromPath("${params.fasta}*", checkIfExists: true).collect()
    : Channel.empty()

    MAKE_METH_BIGWIG(DRAGEN_METHYLATION.out.dragen_output,ch_fasta_reference)
    ch_versions = ch_versions.mix(MAKE_METH_BIGWIG.out.versions)

    emit: 
    dragen_outputs = DRAGEN_METHYLATION.out.dragen_output
    versions = ch_versions

}