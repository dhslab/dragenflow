process TRIM_ADAPTERS {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/dhslab/docker-fastp:251229"

    input:
    tuple val(meta), path(read1file), path(read2file)
    path(adapter1)
    path(adapter2)

    output:
    tuple val(meta), path("*.R1_trimmed.fastq.gz"), path("*.R2_trimmed.fastq.gz"), emit: trimmed_fastqs
    path("versions.yml")   , emit: versions

    script:
    def command_args = [
        read1file  ? "-i ${read1file} -o ${read1file.toString().replaceAll(/\.fq$|\.fastq$|\.fq\.gz$|\.fastq\.gz$/, '')}.R1_trimmed.fastq.gz" : "",
        read2file  ? "-I ${read2file} -O ${read2file.toString().replaceAll(/\.fq$|\.fastq$|\.fq\.gz$|\.fastq\.gz$/, '')}.R2_trimmed.fastq.gz" : "",
        task.cpus  ? "-w ${task.cpus}" : ""
    ].join(' ').trim()
    """
    cat ${adapter1} ${adapter2} > adapters.fa

    fastp -Q -L --adapter_fasta adapters.fa ${command_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        \$(create_fastq_list.py -v)
    END_VERSIONS
    """

    stub:
        def command_args = [
        read1file  ? "-i ${read1file} -o ${read1file.toString().replaceAll(/\.fq$|\.fastq$|\.fq\.gz$|\.fastq\.gz$/, '')}.R1_trimmed.fastq.gz" : "",
        read2file  ? "-I ${read2file} -O ${read2file.toString().replaceAll(/\.fq$|\.fastq$|\.fq\.gz$|\.fastq\.gz$/, '')}.R2_trimmed.fastq.gz" : "",
        task.cpus  ? "-w ${task.cpus}" : ""
    ].join(' ').trim()
    """
    echo fastp ${command_args}
    touch \$(basename ${read1file} .fastq.gz).R1_trimmed.fastq.gz
    touch \$(basename ${read2file} .fastq.gz).R2_trimmed.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 1.0.1
    END_VERSIONS
    """

}
