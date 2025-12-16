include { ANNOTATE_EXPRESSION_TABLES                  } from '../../modules/local/annotate_expression_tables.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW } from '../../modules/local/bedtools_genomecov.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_RV } from '../../modules/local/bedtools_genomecov.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_US } from '../../modules/local/bedtools_genomecov.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD    } from './bedgraph_bedclip_bedgraphtobigwig.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE    } from './bedgraph_bedclip_bedgraphtobigwig.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_UNSTRANDED } from './bedgraph_bedclip_bedgraphtobigwig.nf'

workflow RNASEQ {
    take:
    ch_rnaseq_dragen_output
    ch_rnaseq_transcript_table

    main:
    ch_versions = Channel.empty()

    ANNOTATE_EXPRESSION_TABLES(ch_rnaseq_dragen_output, ch_rnaseq_transcript_table)
    ch_versions = ch_versions.mix(ANNOTATE_EXPRESSION_TABLES.out.versions)
    
    emit:
    dragen_output = DRAGEN_RNA.out.dragen_output
    versions = ch_versions
}

//
// Processes
//

process GET_SIZES_FILE {
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils:240229"

    input:
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    output:
    path('*sizes'), emit: sizes

    script:
    """
    sizes_file=\$(basename inputs/*.fa.fai | sed 's/\\.fa\\.fai\$/.sizes/')
    cut -f 1,2 inputs/*.fa.fai > \$sizes_file
    """
}

process GET_STRANDEDNESS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils:240229"

    input:
    tuple val(meta), path(dragen_output)
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    output:
    tuple val(meta), path(dragen_output), env(strandedness), emit: dragen_output

    script:
    """
    strand_type=\$(head -n 1 *.quant_metrics.csv | cut -d',' -f4)

    # Adjust strandedness based on the result
    if [ "\$strand_type" = "ISR" ]; then
        strandedness="reverse"
    elif [ "\$strand_type" = "ISF" ]; then
        strandedness="forward"
    elif [ "\$strand_type" = "IU" ]; then
        strandedness="unstranded"
    fi
    """
}