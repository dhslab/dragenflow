process CONCATENATE_FASTQLISTS {
    tag "$meta.id"
    label 'process_local'

    input:
    tuple val(meta), path(fastqlists)

    output:
    tuple val(meta), path("${meta.id}.fastq_list.csv"), emit: fastqlist

    script:
    """
    echo "RGID,RGSM,RGLB,Lane,RGPU,Read1File,Read2File" > ${meta.id}.fastq_list.csv
    cat ${fastqlists.join(' ')} >> ${meta.id}.fastq_list.csv
    """
}