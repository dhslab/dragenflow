include { DRAGEN_MULTIALIGN as DRAGEN_FIVE_BASE } from '../../modules/local/dragen_multialign.nf'
include { ANNOTATE_SMALLVARIANTS                } from '../../modules/local/annotate_smallvariants.nf'

workflow FIVE_BASE {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()

    DRAGEN_FIVE_BASE(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_FIVE_BASE.out.versions)

    ch_dragen_output = DRAGEN_FIVE_BASE.out.dragen_output
    ch_fasta = Channel.fromPath(params.fasta).collect()
    ch_vepcache = Channel.fromPath(params.vepcache).collect()

    ANNOTATE_SMALLVARIANTS(ch_dragen_output, ch_fasta, ch_vepcache)
    ch_versions = ch_versions.mix(ANNOTATE_SMALLVARIANTS.out.versions)

    emit:
    versions = ch_versions
}
