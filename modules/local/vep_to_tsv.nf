process VEP_TO_TSV {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-baseimage:latest"

    input:
    tuple val(meta), path(input)
    val(type)

    output:
    tuple val(meta), path("*.tsv"), emit: vep_tsv
    path "versions.yml", emit: versions

    script:
    def args = [
     type.toLowerCase() == "vcf" ? "-v" : "",
     type.toLowerCase() == "cnv" ? "-s" : "",
     type.toLowerCase() == "sv" ? "-s" : "",
     type.toLowerCase() != "cnv" && param.workflow == "somatic" ? "-i 1" : "-i 0"
    ].join(' ').trim()

    def output = input.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')
    """
    vep2table.py $args $input -o $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}