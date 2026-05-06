process VEP_TO_TSV {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-baseimage:latest"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(input)
    val(type)

    output:
    tuple val(meta), path("*.tsv"), emit: vep_tsv
    path "versions.yml", emit: versions

    script:
    def t = type.toLowerCase()
    def isSomatic = params.workflow == "somatic"
    def args = ""
    switch (t) {
        case 'vcf':
            args = "-i ${isSomatic ? '0 1' : '0'} -v"
            break
        case 'sv':
            args = "-i ${isSomatic ? '1' : '0'} -s"
            break
        case 'cnv':
            args = "-i 0 -s"
            break
    }
    def vcf = input.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def output = vcf.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')

    """
    vep2table.py $args $vcf -o $output
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def t = type.toLowerCase()
    def isSomatic = params.workflow == "somatic"
    def args = ""
    switch (t) {
        case 'vcf':
            args = "-i ${isSomatic ? '0 1' : '0'} -v"
            break
        case 'sv':
            args = "-i ${isSomatic ? '1' : '0'} -s"
            break
        case 'cnv':
            args = "-i 0 -s"
            break
    }
    def vcf = input.find{ it ==~ /.*\.(vcf.gz)$/ } ?: ""
    def output = vcf.getName().replaceFirst('\\.vcf\\.gz\$', '.tsv')

    """
    touch "${output}"

    cat <<-END_CMDS > "${task.process}_cmds.txt"
    vep2table.py $args $vcf -o $output
    END_CMDS

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}