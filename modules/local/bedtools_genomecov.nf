process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_high'
    container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_output)
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')
    val(process_type)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"                , emit: versions        

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if( (process_type == 'forward' || process_type == 'reverse') )
        """
        bedtools \\
            genomecov \\
            -ibam *cram \\
            $args \\
            | bedtools sort > ${prefix}.bedGraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """

    else if( process_type == 'unstranded' )
        """
        bedtools \\
            genomecov \\
            -ibam *cram \\
            $args \\
            | bedtools sort > ${meta.id}.unstranded.bedGraph

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
}