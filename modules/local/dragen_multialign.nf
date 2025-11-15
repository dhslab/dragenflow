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
    def input = ""
    if (type == 'fastq') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-fastq-list ${meta.id}.fastq_list.csv --tumor-fastq-list-sample-id ${meta.id}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--fastq-list ${meta.id}.fastq_list.csv --fastq-list-sample-id ${meta.id}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-fastq-list ${meta.tumor}.fastq_list.csv --tumor-fastq-list-sample-id ${meta.tumor} --fastq-list ${meta.normal}.fastq_list.csv --fastq-list-sample-id ${meta.normal}"
        }
    } else if (type == 'cram') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-cram-input ${meta.cram}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--cram-input ${meta.cram}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-cram-input ${meta.tumor} --cram-input ${meta.normal}"
        }
        input += " --cram-reference inputs/${dragen_inputs.input_cram_reference}"
    }
    if (type == 'bam') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-bam-input ${meta.bam}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--bam-input ${meta.bam}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-bam-input ${meta.tumor} --bam-input ${meta.normal}"
        }

    }

    def args_license = task.ext.dragen_license_args ?: ''

    def alignment_params = [
        task.ext.dragen_license_args                  ?: "",
        meta.sex?.toLowerCase() in ['male', 'female'] ? "--sample-sex ${meta.sex}"                                           : "",
        reference_dir                                 ? "--ref-dir ${reference_dir}"                                         : "",
        adapter1                                      ? "--trim-adapter-read1 ${adapter1}"                                   : "",
        adapter2                                      ? "--trim-adapter-read2 ${adapter2}"                                   : "",
        hotspots                                      ? "--vc-somatic-hotspots ${hotspots.min{ it.toString().length() }}"    : "",
        cram_reference                                ? "--cram-reference ${cram_reference.min{ it.toString().length() }}"   : "",
        sv_noise_file                                 ? "--sv-systematic-noise ${sv_noise_file}"                             : "",
        snv_noise_file                                ? "--vc-systematic-noise ${snv_noise_file}"                            : "",
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"               : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"         : "",
        cnv_population_vcf                            ? "--cnv-population-b-allele-vcf ${cnv_population_vcf}"                : "",
        params.dragen_cnv_filter_length               ? "--cnv-filter-length ${params.dragen_cnv_filter_length}"             : "",
        params.dux4caller                             ? "--enable-dux4-caller true --enable-ploidy-estimator true"           : "",
        params.dragen_cnv_merge_distance              ? "--cnv-merge-distance ${params.dragen_cnv_merge_distance}"           : "",
        tandem_duplications                           ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}" : ""
    ].join(' ').trim()

    def dragen_mode_args = [
        intermediate_directory                        ? "--intermediate-results-dir ${intermediate_directory}"               : "",
        intermediate_directory_value                  ? "--intermediate-results-dir ${intermediate_directory_value}"         : "",
        params.alignment_file_format ? "--output-format ${params.alignment_file_format}" : "",
        adapter1 && adapter2 ? "--read-trimmers adapter --trim-adapter-read1 ${adapter1} --trim-adapter-read2 ${adapter2}" : "",
        params.umi ? "--umi-enable true --umi-min-supporting-reads ${params.readfamilysize} --umi-library-type ${params.umi}" : "",
        params.umi && target_bed ? "--umi-metrics-interval-file ${target_bed}" : "",
        "--enable-duplicate-marking ${params.mark_duplicates}",
        params.umi && params.solid_tumor ? "--vc-enable-umi-solid true" : "",
        params.umi && params.liquid_tumor ? "--vc-enable-umi-liquid true" : "",
        hotspots ? "--vc-somatic-hotspots ${hotspots}" : "",
        tandem_duplications ? "--sv-somatic-ins-tandup-hotspot-regions-bed ${tandem_duplications}" : "",
        params.dux4caller ? " --enable-dux4-caller true" : "",
        nirvana_path ? "--enable-variant-annotation true --variant-annotation-assembly ${params.nirvana_assembly} --variant-annotation-data ${nirvana_path}",
        nirvana_path ? "--vc-enable-germline-tagging true" : "--vc-skip-germline-tagging true",
        meta.sex != null ? "--sample-sex ${meta.sex}" : ""
    ].join(' ').trim()
    

    """
    mkdir dragen && \\
    ${task.ext.dragen_exe_path}/dragen -r inputs/${dragen_inputs.reference} ${specified_sex} ${input} ${intermediate_dir} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                
                --output-directory ./dragen --force --output-file-prefix ${meta.id} ${dragen_mode_args} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${task.ext.dragen_exe_path}/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def input = ""
    if (type == 'fastq') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-fastq-list ${meta.id}.fastq_list.csv --tumor-fastq-list-sample-id ${meta.id}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--fastq-list ${meta.id}.fastq_list.csv --fastq-list-sample-id ${meta.id}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-fastq-list ${meta.tumor}.fastq_list.csv --tumor-fastq-list-sample-id ${meta.tumor} --fastq-list ${meta.normal}.fastq_list.csv --fastq-list-sample-id ${meta.normal}"
        }
    } else if (type == 'cram') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-cram-input ${meta.cram}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--cram-input ${meta.cram}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-cram-input ${meta.tumor} --cram-input ${meta.normal}"
        }
        input += " --cram-reference inputs/${dragen_inputs.input_cram_reference}"
    }
    if (type == 'bam') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-bam-input ${meta.bam}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--bam-input ${meta.bam}"
        } else if (params.workflow == "somatic"){
            input = "--tumor-bam-input ${meta.tumor} --bam-input ${meta.normal}"
        }

    }

    def intermediate_dir = task.ext.intermediate_dir ? "--intermediate-results-dir ${task.ext.intermediate_dir}" : ""
    def args_license = task.ext.dragen_license_args ?: ''
    def specified_sex = meta.sex != null ? "--sample-sex ${meta.sex}" : ""
    def tandup_bed = dragen_inputs.tandem_dup_hotspot_bed != null ? "--sv-somatic-ins-tandup-hotspot-regions-bed inputs/${dragen_inputs.tandem_dup_hotspot_bed}" : ""
    def dux4caller = params.dux4caller == true ? " --enable-dux4-caller true" : ""
    def hotspotvcf = dragen_inputs.hotspot_vcf != null ? "--vc-somatic-hotspots inputs/${dragen_inputs.hotspot_vcf}" : ""

    def dragen_mode_args = ""
    def args = params.dragen_args ?: ""

    if (params.udiumi == true && dragen_inputs.target_bed_file){
        dragen_mode_args = "--umi-enable true --umi-min-supporting-reads ${params.readfamilysize} --umi-library-type random-simplex --umi-metrics-interval-file inputs/${dragen_inputs.target_bed_file}"
        
        if (params.solid_tumor){
            dragen_mode_args += " --vc-enable-umi-solid true"
        } else {
            dragen_mode_args += " --vc-enable-umi-liquid true"
        }

    } else {
        dragen_mode_args = "--enable-duplicate-marking ${params.mark_duplicates} --read-trimmers adapter --trim-adapter-read1 inputs/${dragen_inputs.dragen_adapter1} --trim-adapter-read2 inputs/${dragen_inputs.dragen_adapter2}"
    }

    if (params.workflow == "rna"){
        def downsampleargs = params.downsample_rna ? " --enable-down-sampler true --down-sampler-reads 200000000" : ""
        dragen_mode_args += " --enable-variant-caller true --enable-rna true -a inputs/${dragen_inputs.annotation_file} --rrna-filter-enable true --enable-rna-quantification true --enable-rna-gene-fusion true ${downsampleargs}"
    
    } else if (params.workflow == "5mc"){
        dragen_mode_args += " --enable-methylation-calling true --methylation-protocol directional --methylation-generate-cytosine-report true --methylation-compress-cx-report true"
        
    } else if (params.workflow == "tumor" || params.workflow == "somatic"){
        dragen_mode_args += " --enable-variant-caller true --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 3 --sv-min-candidate-spanning-count 3 --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} --sv-enable-somatic-ins-tandup-hotspot-regions true ${tandup_bed}"
        if (params.target_bed_file){
            dragen_mode_args += "  --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-population-b-allele-vcf inputs/${dragen_inputs.pop_af_vcf} --cnv-enable-self-normalization true --cnv-target-bed inputs/${dragen_inputs.target_bed_file} --sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += " --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-population-b-allele-vcf inputs/${dragen_inputs.pop_af_vcf}${dux4caller}"
        }

        if (params.nirvana){
            dragen_mode_args += " --vc-enable-germline-tagging true"
        } else {
            dragen_mode_args += " --vc-skip-germline-tagging true"
        }

    } else if (params.workflow == "germline"){
        dragen_mode_args += " --enable-variant-caller true --enable-sv true --sv-output-contigs true --sv-use-overlap-pair-evidence true"
        if (params.target_bed_file){
            dragen_mode_args += " --sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += " --enable-cnv true --cnv-enable-self-normalization true"
        }

    } else if (params.workflow == "align"){
        dragen_mode_args += " --enable-variant-caller false --enable-sv false --enable-cnv false"
    }

    if (params.nirvana){
        dragen_mode_args += " --enable-variant-annotation true --variant-annotation-assembly GRCh38 --variant-annotation-data inputs/${dragen_inputs.nirvana}"
    }

    """
    mkdir dragen && \\
    echo ${task.ext.dragen_exe_path}/dragen -r inputs/${dragen_inputs.reference} ${specified_sex} ${input} ${intermediate_dir} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps true \\
                --gc-metrics-enable true \\
                --output-format ${params.alignment_file_format} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} ${dragen_mode_args} ${args}
    cp ${params.dragendir}/* dragen/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}
