process MAKE_METH_BIGWIG {
    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_output)
    path(reference)

    output:
    path("*.meth.bw"), emit: bigwig
    path "versions.yml"    , emit: versions        

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chromsizes = reference.find{ it ==~ /.*\.(fai)$/ }?.with { "$it" } 

    """
    dragen_meth_to_bw.py -s ${chromsizes} -o ${meta.id}.meth.bw ${meta.id}.CX_report.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_meth_to_bw.py: \$(dragen_meth_to_bw.py --version | cut -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chromsizes = reference.find{ it ==~ /.*\.(fai)$/ }?.with { "$it" } 

    """
    touch ${meta.id}.meth.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_meth_to_bed.py: \$(dragen_meth_to_bed.py --version | cut -f 2)
    END_VERSIONS
    """

}