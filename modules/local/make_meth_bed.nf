process MAKE_METH_BED {
    tag "$meta.id"
    label 'process_medium'
    container "ghcr.io/dhslab/baseimage:250222"

    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_output)

    output:
    path("*.meth.bed.gz"), emit: bed
    path "versions.yml"        , emit: versions        

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    dragen_meth_to_bed.py -o ${meta.id}.meth.bed.gz ${meta.id}.CX_report.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_meth_to_bed.py: \$(dragen_meth_to_bed.py --version | cut -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${meta.id}.meth.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_meth_to_bed.py: \$(dragen_meth_to_bed.py --version | cut -f 2)
    END_VERSIONS
    """
}