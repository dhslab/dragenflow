process MAKE_FASTQLIST {
    tag "$meta.id"
    label 'process_tiny'
    container "ghcr.io/dhslab/docker-python3:231224"

    input:
    tuple val(meta), path(inputs)
    val(type)

    output:
    tuple val(meta), path("${meta.id}.fastq_list.csv"), emit: fastqlist
    path 'versions.yml', emit: versions

    script:
    def input_flags = ""

    if (type=="demux_path"){
        input_flags = "-d ${inputs}"

    } else if (type=="fastq_list"){
        input_flags = "-f ${inputs}"
   
    } else if (type=="fastq" || type=="reads" || type=="fastqs"){
        input_flags = "-1 ${inputs[0]} -2 ${inputs[1]}"

    } else {
        error "Unknown input type: ${type}"
    }

    """
    make_fastqlist.py -i ${meta.id} ${input_flags}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(make_fastqlist.py -v)
    END_VERSIONS
    """

    stub:
    def input_flags = ""
    if (type=="demux_path"){
        input_flags = "-d ${inputs}"

    } else if (type=="fastq_list"){
        input_flags = "-f ${inputs}"
   
    } else if (type=="fastq" || type=="reads"){
        input_flags = "-1 ${inputs[0]} -2 ${inputs[1]}"

    } else {
        error "Unknown input type: ${type}"
    }

    """
    make_fastqlist.py -i ${meta.id} ${input_flags}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(make_fastqlist.py -v)
        container: ${task.container}
    END_VERSIONS
    """

}
