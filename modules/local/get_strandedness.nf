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