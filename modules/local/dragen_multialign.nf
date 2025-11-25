process DRAGEN_MULTIALIGN {
    tag "${meta.id}"
    label 'dragen'
    container "${task.ext.dragen_container}"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename.split('/')[1] }, mode:'copy'

//    input:
//    tuple val(meta), val(type), path("*")
//    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

    input:
    tuple val(meta), path(reads, stageAs: "fastq_files/*"), path(fastq_list), path(alignment_file)
    tuple val(intermediate_directory_value), path(intermediate_directory)
    path(reference_dir)
    path(dbsnp)
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
    path("dragen/${meta.id}_usage.txt"), emit: usage, optional: true
    path "versions.yml",    emit: versions

    when:
    params.run_dragen == true

    script:
    def exe_path = "${task.ext.dragen_path}"
    def input = ""
    if (params.workflow == "rna" || params.workflow == "tumor"){
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--tumor-bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--tumor-cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "bsseq" || params.workflow == "germline" || params.workflow == "align") {
        input = [
            alignment_file.find{ it ==~ /.*\.bam$/  }?.with{ "--bam-input ${it}"  }                                        ?:
            alignment_file.find{ it ==~ /.*\.cram$/ }?.with{ "--cram-input ${it}" }                                        ?:
            fastq_list.toString().endsWith('csv')    ? "--fastq-list ${fastq_list} --fastq-list-sample-id ${meta.id}" :
            error("Input file is not a BAM, CRAM, or CSV file.")
        ].join(' ').trim()

    } else if (params.workflow == "somatic"){
        input = [
            fastq_list.toString().endsWith('csv')    ? "--tumor-fastq-list ${fastq_list} --tumor-fastq-list-sample-id ${meta.tumor_id} --fastq-list ${fastq_list} --fastq-list-sample-id ${meta.normal_id}" :
            error("Input file is not a CSV file.")
        ].join(' ').trim()
    }

    def alignment_params = [
        task.ext.dragen_args                          ?: "",
        params.extra_dragen_args                      ?: "",
        task.ext.dragen_license_args                  ?: "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"    : "",
        dbsnp                                         ? "--dbsnp ${dbsnp}"                                                    : "",
        params.alignment_file_format                  ? "--output-format ${params.alignment_file_format}"                     : "",
        meta.sex?.toLowerCase() in ['male', 'female'] ? "--sample-sex ${meta.sex}"                                            : "",
        adapter1 && adapter2                          ? "--read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2}" : "",
        params.umi                                    ? "--umi-enable true --umi-library-type=${params.umi}"                  : "",
        params.umi && params.readfamilysize           ? "--umi-min-supporting-reads ${readfamilysize}"                        : "",
        params.umi && target_bed                      ? "--umi-metrics-interval-file ${target_bed}"                           : "",
        params.umi && liquid_tumor                    ? "--vc-enable-umi-liquid true"                                         : "",
        params.umi && solid_tumor                     ? "--vc-enable-umi-solid true"                                          : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"     : "",
        params.variant_caller && target_bed           ? "--vc-target-bed ${target_bed}"                                       : "",              
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        params.sv_caller && target_bed                ? "--sv-target-bed ${target_bed} --sv-exome true"                       : "",
        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                 : "",
        params.cnv_caller && params.dragen_cnv_filter_length               ? "--cnv-filter-length ${params.dragen_cnv_filter_length}"              : "",
        params.cnv_caller && params.dragen_cnv_merge_distance              ? "--cnv-merge-distance ${params.dragen_cnv_merge_distance}"            : "",
        tandem_duplications                           ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}"  : "",
        params.use_nirvana && nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        params.use_nirvana && nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
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

    # Copy and rename DRAGEN usage
    find dragen/ \\
        -type f \\
        -name "*_usage.txt" \\
        -exec mv "{}" "dragen/${meta.id}_usage.txt" \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${task.ext.dragen_exe_path}/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = "${task.ext.dragen_path}"
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
        task.ext.dragen_args                          ?: "",
        task.ext.dragen_license_args                  ?: "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"                : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"          : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                          : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"    : "",
        dbsnp                                         ? "--dbsnp ${dbsnp}"                                                    : "",
        params.alignment_file_format                  ? "--output-format ${params.alignment_file_format}"                     : "",
        meta.sex?.toLowerCase() in ['male', 'female'] ? "--sample-sex ${meta.sex}"                                            : "",
        adapter1 && adapter2                          ? "--read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2}" : "",
        params.umi                                    ? "--umi-enable true --umi-library-type=${params.umi}"                  : "",
        params.umi && params.readfamilysize           ? "--umi-min-supporting-reads ${readfamilysize}"                        : "",
        params.umi && target_bed                      ? "--umi-metrics-interval-file ${target_bed}"                           : "",
        params.umi && liquid_tumor                    ? "--vc-enable-umi-liquid true"                                         : "",
        params.umi && solid_tumor                     ? "--vc-enable-umi-solid true"                                          : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"     : "",
        params.variant_caller && target_bed           ? "--vc-target-bed ${target_bed}"                                       : "",              
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                             : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                              : "",
        params.sv_caller && target_bed                ? "--sv-target-bed ${target_bed} --sv-exome true"                       : "",
        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                 : "",
        params.cnv_caller && params.dragen_cnv_filter_length               ? "--cnv-filter-length ${params.dragen_cnv_filter_length}"              : "",
        params.cnv_caller && params.dragen_cnv_merge_distance              ? "--cnv-merge-distance ${params.dragen_cnv_merge_distance}"            : "",
        tandem_duplications                           ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}"  : "",
        params.use_nirvana && nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}" : "",
        params.use_nirvana && nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true"
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
    cp ${projectDir}/assets/stub/dragen_path/* dragen/
    touch dragen/${meta.id}_usage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}