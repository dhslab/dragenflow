include { DRAGEN_MULTIALIGN as DRAGEN_TUMOR     } from '../../modules/local/dragen_multialign.nf'
include { ANNOTATE_SMALLVARIANTS                } from '../../modules/local/annotate_smallvariants.nf'

workflow TUMOR {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_TUMOR(input_data, dragen_inputs)
    ch_dragen_output = ch_dragen_output.mix(DRAGEN_TUMOR.out.dragen_output)
    ch_versions = ch_versions.mix(DRAGEN_TUMOR.out.versions)

    ANNOTATE_SMALLVARIANTS(ch_dragen_output, Channel.fromPath(params.fasta), Channel.fromPath(params.vepcache))
    ch_versions = ch_versions.mix(ANNOTATE_SMALLVARIANTS.out.versions)

    emit:
    dragen_output = ch_dragen_output
    ch_versions
}
