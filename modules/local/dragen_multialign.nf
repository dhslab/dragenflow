process DRAGEN_MULTIALIGN {
    tag "${meta.id}"
    label 'dragen'
    label 'dragenalign'
    container "${ext.dragen_aws_image}" ?: "${params.dragen_container}"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename.split('/')[1] }, mode:'copy'

//    input:
//    tuple val(meta), val(type), path("*")
//    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    input:
    tuple val(meta), path(reads, stageAs: "fastq_files/*"), path(fastq_list), path(alignment_file)
    tuple val(intermediate_directory_value), path(intermediate_directory)
    path(reference_dir)
    path(adapter1)
    path(adapter2)
    path(cram_reference)
    path(sv_noise_file)
    path(snv_noise_file)
    path(hotspots)
    path(cnv_population_vcf)
    path(tandem_duplications)
    path(target_bed)
    path(annotation_gtf)
    path(nirvana_path)

    output:
    tuple val(meta), path ("dragen/*"), emit: dragen_output
    path "versions.yml",    emit: versions

    script:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "${params.aws_dragen_path}" : "${params.local_dragen_path}"
    def input = ""
    if (params.workflow == "rna" || params.workflow == "tumor"){
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--tumor-bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--tumor-cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--fastq-list ${fastq_list} --fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "somatic"){
        input = [
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.tumor_id}" :
            fastq_list.toString().endsWith('csv')    ? "--fastq-list ${fastq_list} --fastq-list-sample-id ${meta.normal_id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()
    }

    def alignment_params = [
        task.ext.dragen_license_args                  ?: "",
        meta.sex?.toLowerCase() in ['male', 'female'] ? "--sample-sex ${meta.sex}"                                            : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        adapter1                                      ? "--trim-adapter-read1 ${adapter1}"                                    : "",
        adapter2                                      ? "--trim-adapter-read2 ${adapter2}"                                    : "",
        params.alignment_file_format ? "--output-format ${params.alignment_file_format}"                                      : "",
        params.umi ? "--umi-enable true --umi-min-supporting-reads ${params.readfamilysize} --umi-library-type ${params.umi}" : "",
        params.umi && target_bed ? "--umi-metrics-interval-file ${target_bed}"                                                : "",
        params.mark_duplicates ? "--enable-duplicate-marking true"                                                            : "--enable-duplicate-marking false",
        params.umi && params.solid_tumor ? "--vc-enable-umi-solid true"                                                       : "",
        params.umi && params.liquid_tumor ? "--vc-enable-umi-liquid true"                                                     : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"     : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"    : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        params.dragen_cnv_filter_length               ? "--cnv-filter-length ${params.dragen_cnv_filter_length}"              : "",
        params.dux4caller                             ? "--enable-dux4-caller true --enable-ploidy-estimator true"            : "",
        params.dragen_cnv_merge_distance              ? "--cnv-merge-distance ${params.dragen_cnv_merge_distance}"            : "",
        tandem_duplications                           ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}"  : "",
        nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
    ].join(' ').trim()

    """
    mkdir dragen && \\
    ${exe_path}/bin/dragen \\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                ${alignment_params} \\
                ${input} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${task.ext.dragen_exe_path}/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = ['dragenaws', 'awsbatch'].any{ workflow.profile.contains(it) } ? "${params.aws_dragen_path}" : "${params.local_dragen_path}"
    def input = ""
    if (params.workflow == "rna" || params.workflow == "tumor"){
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--tumor-bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--tumor-cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--fastq-list ${fastq_list} --fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "somatic"){
        input = [
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.tumor_id}" :
            fastq_list.toString().endsWith('csv')    ? "--fastq-list ${fastq_list} --fastq-list-sample-id ${meta.normal_id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()
    }

    def alignment_params = [
        task.ext.dragen_license_args                  ?: "",
        meta.sex?.toLowerCase() in ['male', 'female'] ? "--sample-sex ${meta.sex}"                                            : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        adapter1                                      ? "--trim-adapter-read1 ${adapter1}"                                    : "",
        adapter2                                      ? "--trim-adapter-read2 ${adapter2}"                                    : "",
        params.alignment_file_format ? "--output-format ${params.alignment_file_format}"                                      : "",
        params.umi ? "--umi-enable true --umi-min-supporting-reads ${params.readfamilysize} --umi-library-type ${params.umi}" : "",
        params.umi && target_bed ? "--umi-metrics-interval-file ${target_bed}"                                                : "",
        params.mark_duplicates ? "--enable-duplicate-marking true"                                                            : "--enable-duplicate-marking false",
        params.umi && params.solid_tumor ? "--vc-enable-umi-solid true"                                                       : "",
        params.umi && params.liquid_tumor ? "--vc-enable-umi-liquid true"                                                     : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"     : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"    : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        params.dragen_cnv_filter_length               ? "--cnv-filter-length ${params.dragen_cnv_filter_length}"              : "",
        params.dux4caller                             ? "--enable-dux4-caller true --enable-ploidy-estimator true"            : "",
        params.dragen_cnv_merge_distance              ? "--cnv-merge-distance ${params.dragen_cnv_merge_distance}"            : "",
        tandem_duplications                           ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}"  : "",
        nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
    ].join(' ').trim()

    """
    mkdir dragen && \\
    echo ${exe_path}/bin/dragen \\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                ${alignment_params} \\
                ${input} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} > dragen/command.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}