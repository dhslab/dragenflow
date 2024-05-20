process DRAGEN_MULTIALIGN {
    tag "${meta.id}"
    label 'dragen'
    label 'dragenalign'
    container "${ext.dragen_aws_image}" ?: "${params.dragen_container}"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename == "versions.yml" ? null : filename.split('/')[1] }, mode:'copy'

    input:
    tuple val(meta), val(type), path("*")
    tuple val(dragen_inputs), path("*", stageAs: 'inputs/*')

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
            input = "--tumor-fastq-list ${meta.id}.fastq_list.csv --tumor-fastq-list-sample-id ${meta.tumor} --fastq-list ${meta.id}.fastq_list.csv --fastq-list-sample-id ${meta.normal}"
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

    def dragen_mode_args = ""

    if (params.workflow == "rna"){
        def downsampleargs = params.downsample_rna ? " --enable-down-sampler true --down-sampler-reads 100000000" : ""
        dragen_mode_args = "--enable-variant-caller true --enable-rna true -a inputs/${dragen_inputs.annotation_file} --rrna-filter-enable true --enable-rna-quantification true --enable-rna-gene-fusion true ${downsampleargs}"
    
    } else if (params.workflow == "5mc"){
        dragen_mode_args = "--enable-methylation-calling true --methylation-protocol directional --methylation-generate-cytosine-report true --methylation-compress-cx-report true"
        
    } else if (params.workflow == "tumor" || params.workflow == "somatic"){
        def tandup_bed = dragen_inputs.tandem_dup_hotspot_bed != null ? "--sv-somatic-ins-tandup-hotspot-regions-bed inputs/${dragen_inputs.tandem_dup_hotspot_bed}" : ""
        def dux4caller = params.dux4caller == true ? " --enable-dux4-caller true" : ""
        def hotspotvcf = dragen_inputs.hotspot_vcf != null ? "--vc-somatic-hotspots inputs/${dragen_inputs.hotspot_vcf}" : ""
        dragen_mode_args = "--enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} ${hotspotvcf} --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} --sv-enable-somatic-ins-tandup-hotspot-regions true ${tandup_bed}"
        if (params.targeted_sequencing == true || params.target_bed_file){
            dragen_mode_args += " --sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += " --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-population-b-allele-vcf inputs/${dragen_inputs.pop_af_vcf}${dux4caller}"
        }

    } else if (params.workflow == "germline"){
        dragen_mode_args = "--enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} --enable-sv true --sv-output-contigs true --sv-use-overlap-pair-evidence true"
        if (params.targeted_sequencing == true || params.target_bed_file){
            dragen_mode_args += "--sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += "--enable-cnv true --cnv-enable-self-normalization true --enable-cyp2b6 true --enable-cyp2d6 true --enable-gba true --enable-smn true --repeat-genotype-enable true"
        }
    } else if (params.workflow == "align"){
        dragen_mode_args = "--enable-variant-caller false --enable-sv false"
    }

    """
    mkdir dragen && \\
    /opt/edico/bin/dragen -r inputs/${dragen_inputs.reference} ${specified_sex} ${input} ${intermediate_dir} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps=true \\
                --gc-metrics-enable true \\
                --enable-duplicate-marking ${params.mark_duplicates} \\
                --read-trimmers adapter \\
                --trim-adapter-read1 inputs/${dragen_inputs.dragen_adapter1} \\
                --trim-adapter-read2 inputs/${dragen_inputs.dragen_adapter2} \\
                --output-format ${params.alignment_file_format} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} ${dragen_mode_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def input = ""
    if (type == 'fastq') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-fastq-list fastq_list.csv --tumor-fastq-list-sample-id ${meta.id}"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--fastq-list fastq_list.csv --fastq-list-sample-id ${meta.id}"
        }
    } else if (type == 'cram') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-cram-input *.cram"
        } else if (params.workflow == "5mc" || params.workflow == "germline") {
            input = "--cram-input *.cram"
        }
    }
    if (type == 'bam') {
        if (params.workflow == "rna" || params.workflow == "tumor"){
            input = "--tumor-bam-input *.bam"
        } else if (params.workflow == "5mc" || params.workflow == "germline" || params.workflow == "align") {
            input = "--bam-input *.bam"
        }
    }
    def intermediate_dir = task.ext.intermediate_dir ? "--intermediate-results-dir ${task.ext.intermediate_dir}" : ""
    def args_license = task.ext.dragen_license_args ?: ''
    def specified_sex = meta.sex != null ? "--sample-sex ${meta.sex}" : ""

    def dragen_mode_args = ""

    if (params.workflow == "rna"){
        def downsampleargs = params.downsample_rna ? " --enable-down-sampler true --down-sampler-reads 100000000" : ""
        dragen_mode_args = "--enable-variant-caller true --enable-rna true -a inputs/${dragen_inputs.annotation_file} --rrna-filter-enable true --enable-rna-quantification true --enable-rna-gene-fusion true ${downsampleargs}"
    
    } else if (params.workflow == "5mc"){
        dragen_mode_args = "--enable-methylation-calling true --methylation-protocol directional --methylation-generate-cytosine-report true"
        
    } else if (params.workflow == "tumor"){
        def tandup_bed = dragen_inputs.tandem_dup_hotspot_bed != null ? "--sv-somatic-ins-tandup-hotspot-regions-bed inputs/${dragen_inputs.tandem_dup_hotspot_bed}" : ""
        def dux4caller = params.dux4caller == true ? " --enable-dux4-caller true" : ""
        def hotspotvcf = params.hotspot_vcf != null ? "--vc-somatic-hotspots inputs/${dragen_inputs.hotspot_vcf}" : ""
        dragen_mode_args = "--enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} ${hotspotvcf} --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} --sv-enable-somatic-ins-tandup-hotspot-regions true ${tandup_bed}"

        dragen_mode_args = "--enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} ${hotspotvcf} --vc-systematic-noise inputs/${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 --sv-use-overlap-pair-evidence true --sv-systematic-noise inputs/${dragen_inputs.sv_noisefile} --sv-enable-somatic-ins-tandup-hotspot-regions true ${tandup_bed}"
        if (params.targeted_sequencing == true || params.target_bed_file){
            dragen_mode_args += " --sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += " --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-population-b-allele-vcf inputs/${dragen_inputs.pop_af_vcf}${dux4caller}"
        }

    } else if (params.workflow == "germline"){
        dragen_mode_args = "--enable-variant-caller true --dbsnp inputs/${dragen_inputs.dbsnp} --enable-sv true --sv-output-contigs true --sv-use-overlap-pair-evidence true"
        if (params.targeted_sequencing == true || params.target_bed_file){
            dragen_mode_args += "--sv-exome true --sv-call-regions-bed inputs/${dragen_inputs.target_bed_file} --vc-target-bed inputs/${dragen_inputs.target_bed_file}"
        } else {
            dragen_mode_args += "--enable-cnv true --cnv-enable-self-normalization true --enable-cyp2b6 true --enable-cyp2d6 true --enable-gba true --enable-smn true --repeat-genotype-enable true"
        }
    } else if (params.workflow == "align"){
        dragen_mode_args = "--enable-variant-caller false --enable-sv false"
    }

    """
    mkdir dragen && \\
    echo /opt/edico/bin/dragen -r inputs/${dragen_inputs.reference} ${input} ${intermediate_dir} ${args_license}\\
                --enable-map-align true \\
                --enable-sort true \\
                --enable-bam-indexing true \\
                --enable-map-align-output true \\
                --qc-coverage-ignore-overlaps=true \\
                --gc-metrics-enable true \\
                --enable-duplicate-marking ${params.mark_duplicates} \\
                --read-trimmers adapter \\
                --trim-adapter-read1 inputs/${dragen_inputs.dragen_adapter1} \\
                --trim-adapter-read2 inputs/${dragen_inputs.dragen_adapter2} \\
                --output-format ${params.alignment_file_format} \\
                --output-directory ./dragen --force --output-file-prefix ${meta.id} ${dragen_mode_args} > ./dragen/${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat $projectDir/assets/stub/versions/dragen_version.txt)
    END_VERSIONS
    """
}
