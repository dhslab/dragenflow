process ANNOTATE_RNASEQ {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-python3:240301"

    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_output)
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    output:
    tuple val(meta), path("${meta.id}.quant.genes.annotated.tsv"), emit: genes
    tuple val(meta), path("${meta.id}.quant.annotated.tsv"), emit: transcripts
    path "versions.yml",    emit: versions

    script:
    """
    add_genename_todragenrna.py -i ${meta.id}.quant.genes.sf -d inputs/${dragen_inputs.transcript_table} -o ${meta.id}.quant.genes.annotated.tsv && \\
    add_genename_todragenrna.py -t -i ${meta.id}.quant.sf -d inputs/${dragen_inputs.transcript_table} -o ${meta.id}.quant.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_genename_todragenrna.py: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.quant.genes.annotated.tsv"
    touch "${meta.id}.quant.annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(cat $projectDir/assets/stub/versions/vep_version.yaml)
    END_VERSIONS
    """
}
