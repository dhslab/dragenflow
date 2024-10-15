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