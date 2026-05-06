process CONCATENATE_FASTQLISTS {
    tag "$meta.id"
    label 'process_local'

    input:
    tuple val(meta), path(fastqlists)

    output:
    tuple val(meta), path("${meta.id}.fastq_list.csv"), path("${meta.id}.fastqs.csv"), emit: fastqs

    script:
    """ 
    concatenate_fastqlists.py ${meta.id} ${fastqlists.join(' ')}
    """
}