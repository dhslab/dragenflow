/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS/FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText   } from '../subworkflows/local/utils_nfcore_dragenflow_pipeline'

include { PARSE_INPUT_SAMPLESHEET     } from '../modules/local/parse_input_samplesheet'
include { SAMPLESHEET_CHECK           } from '../modules/local/samplesheet_check.nf'

include { GATHER_ALIGNMENT_SAMPLES      } from '../subworkflows/local/gather_alignment_samples.nf'
include { PREPARE_SOMATIC_FASTQS        } from '../subworkflows/local/utils_somatic_fastq/'
include { MAKE_HOTSPOT_VCF              } from '../modules/local/make_hotspot_vcf.nf'
include { DRAGEN_MULTIALIGN             } from '../modules/local/dragen_multialign.nf'
include { ANNOTATE_VARIANTS             } from '../modules/local/annotate_variants.nf'
include { VEP_TO_TSV as VARIANTS_TO_TSV } from '../modules/local/vep_to_tsv.nf'
include { ANNOTATE_SV_VARIANTS          } from '../modules/local/annotate_sv_variants.nf'
include { VEP_TO_TSV as SV_TO_TSV       } from '../modules/local/vep_to_tsv.nf'
include { ANNOTATE_CNV_VARIANTS         } from '../modules/local/annotate_cnv_variants.nf'
include { VEP_TO_TSV as CNV_TO_TSV      } from '../modules/local/vep_to_tsv'
include { MAKE_METH_BED                 } from '../modules/local/make_meth_bed.nf'
include { MAKE_METH_BIGWIG              } from '../modules/local/make_meth_bigwig.nf'
include { ANNOTATE_EXPRESSION_TABLES    } from '../modules/local/annotate_expression_tables.nf'

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// DRAGEN reference directory
ch_reference_dir = params.refdir
    ? Channel.fromPath(params.refdir, type: 'dir', checkIfExists: true).collect()
    : Channel.empty()

ch_dbsnp = params.dbsnp
    ? Channel.fromPath(params.dbsnp, checkIfExists: true).collect()
    : []

// DRAGEN adapter sequences for read 1
ch_adapter1_file = params.adapter1
    ? Channel.fromPath(params.adapter1, checkIfExists: true).collect()
    : []

// DRAGEN adapter sequences for read 2
ch_adapter2_file = params.adapter2
    ? Channel.fromPath(params.adapter2, checkIfExists: true).collect()
    : []

// DRAGEN intermediate directory
if (params.intermediate_dir?.toString()?.startsWith('/staging')) {
    ch_intermediate_dir = Channel.of(params.intermediate_dir).map{ [ it, [] ] }.collect()
} else if (params.intermediate_dir) {
    ch_intermediate_dir = Channel.fromPath(params.intermediate_dir).map{ [ [], it ] }.collect()
} else {
    ch_intermediate_dir = [ [], [] ]
}

// DRAGEN hotspots
ch_hotspot_vcf = params.hotspot_vcf
    ? Channel.fromPath("${params.hotspot_vcf}*", checkIfExists: true).collect()
    : []

ch_hotspot_bed = params.hotspot_bed
    ? Channel.fromPath(params.hotspot_bed, checkIfExists: true).collect()
    : []

// DRAGEN tandem duplications
ch_dragen_tandem_dup_hotspots = params.dragen_tandem_dup_hotspots
    ? Channel.fromPath(params.dragen_tandem_dup_hotspots, checkIfExists: true).collect()
    : []

// SNV systematic noise BED file
ch_snv_noisefile = params.snv_noisefile
    ? Channel.fromPath(params.snv_noisefile, checkIfExists: true).collect()
    : []

// SV systematic noise BED file
ch_sv_noisefile = params.sv_noisefile
    ? Channel.fromPath(params.sv_noisefile, checkIfExists: true).collect()
    : []

// High confidence CNV VCF file
ch_cnv_population_vcf = params.cnv_population_vcf
    ? Channel.fromPath(params.cnv_population_vcf, checkIfExists: true).collect()
    : []

// CRAM reference file
ch_cram_reference = params.cram_reference
    ? Channel.fromPath("${params.cram_reference}*", checkIfExists: true).collect()
    : []

// Target bed file
ch_target_bed = params.target_bedfile
    ? Channel.fromPath(params.target_bedfile, checkIfExists: true).collect()
    : []

// FastA reference
ch_fasta_reference = params.fasta
    ? Channel.fromPath("${params.fasta}*", checkIfExists: true).collect()
    : Channel.empty()

// Vep cache
ch_vep_cache = params.vep_cache
    ? Channel.fromPath(params.vep_cache, type: 'dir', checkIfExists: true).collect()
    : Channel.empty()

// Cytobands
ch_cytobands = params.cytobands
    ? Channel.fromPath("${params.cytobands}*", checkIfExists: true).collect()
    : []

ch_annotation_gtf = params.annotation_gtf
    ? Channel.fromPath(params.annotation_gtf, checkIfExists: true).collect()
    : []

ch_nirvana_path = params.nirvana_path && params.use_nirvana == true
    ? Channel.fromPath(params.nirvana_path, type: 'dir', checkIfExists: true).collect()
    : []

/*
~~~~~~~~~~~~~~~~~~
MultiQC parameters
~~~~~~~~~~~~~~~~~~
*/

// Config
ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)

// Custom config
ch_multiqc_custom_config = params.multiqc_config
    ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
    : Channel.empty()

// Logo
ch_multiqc_logo = params.multiqc_logo
    ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
    : Channel.empty()

// Methods description
ch_multiqc_custom_methods_description = params.multiqc_methods_description
    ? file(params.multiqc_methods_description, checkIfExists: true)
    : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def generateMetaFromCsv(csv_file) {
    def lines = csv_file.text.readLines()
    def headers = lines[0].split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
    return lines.drop(1).collect { line ->
        def fields = line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/)*.replaceAll(/^"|"$/, '')
        [headers, fields].transpose().collectEntries { k, v -> v ? [(k): v] : [:] }
    }.findAll { it }
}

workflow DRAGENFLOW {

    take:
    ch_input_samplesheet  // channel: [ path(file) ]

    main:
    ch_versions          = Channel.empty()
    ch_dragen_usage      = Channel.empty()
    ch_demux_output      = Channel.empty()
    ch_alignment_samples = Channel.empty() // channel: [ val(meta), path(reads), path(fastq_list), path(alignment_files) ]
    
    //
    // MODULE: Parse input samplesheet to get sample meta for demux, align, and analysis samples. 
    //         Output of this process is 3 csv files, one for each entry point.
    //
    PARSE_INPUT_SAMPLESHEET (
        ch_input_samplesheet
    )
    ch_versions = ch_versions.mix(PARSE_INPUT_SAMPLESHEET.out.versions)

    GATHER_ALIGNMENT_SAMPLES (
        PARSE_INPUT_SAMPLESHEET.out.samples_to_align,
        ch_demux_output.ifEmpty([]),
        ch_cram_reference
    )
    ch_versions          = ch_versions.mix(GATHER_ALIGNMENT_SAMPLES.out.versions)
    ch_alignment_samples = ch_alignment_samples.mix(GATHER_ALIGNMENT_SAMPLES.out.samples)

    if (params.hotspot_bed){
        MAKE_HOTSPOT_VCF(
            ch_hotspot_bed,
            ch_fasta_reference
        )
        ch_versions = ch_versions.mix(MAKE_HOTSPOT_VCF.out.versions)
        ch_hotspot_vcf = MAKE_HOTSPOT_VCF.out.hotspot_vcf
    }

    if (params.workflow == 'somatic') {
        // for somatic workflow, need to assemble tumor and normal data
        PREPARE_SOMATIC_FASTQS(ch_alignment_samples)
        ch_alignment_samples = PREPARE_SOMATIC_FASTQS.out.samples
    }

    ch_alignment_samples.dump(tag:'alignment_samples',pretty:true)

    DRAGEN_MULTIALIGN (
        ch_alignment_samples,
        ch_intermediate_dir,
        ch_reference_dir,
        ch_dbsnp,
        ch_adapter1_file,
        ch_adapter2_file,
        ch_cram_reference,
        ch_sv_noisefile,
        ch_snv_noisefile,
        ch_hotspot_vcf,
        ch_cnv_population_vcf,
        ch_dragen_tandem_dup_hotspots,
        ch_target_bed,
        ch_annotation_gtf,
        ch_nirvana_path
    )
    ch_versions     = ch_versions.mix(DRAGEN_MULTIALIGN.out.versions)
    ch_dragen_usage = ch_dragen_usage.mix(DRAGEN_MULTIALIGN.out.usage)

    if (params.variant_caller == true) {
        ANNOTATE_VARIANTS (
            DRAGEN_MULTIALIGN.out.dragen_output,
            ch_fasta_reference,
            ch_vep_cache
        )
        ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

        VARIANTS_TO_TSV (ANNOTATE_VARIANTS.out.vcf,Channel.value('vcf'))
        ch_versions = ch_versions.mix(VARIANTS_TO_TSV.out.versions)

    }

    if (params.sv_caller == true) {
        ANNOTATE_SV_VARIANTS (
            DRAGEN_MULTIALIGN.out.dragen_output,
            ch_fasta_reference,
            ch_vep_cache,
            ch_cytobands
        )
        ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

        SV_TO_TSV (ANNOTATE_SV_VARIANTS.out.vcf,Channel.value('sv'))
        ch_versions = ch_versions.mix(SV_TO_TSV.out.versions)
    }

    if (params.cnv_caller == true) {
        ANNOTATE_CNV_VARIANTS (
            DRAGEN_MULTIALIGN.out.dragen_output,
            ch_fasta_reference,
            ch_vep_cache,
            ch_cytobands
        )
        ch_versions = ch_versions.mix(ANNOTATE_CNV_VARIANTS.out.versions)

        CNV_TO_TSV (ANNOTATE_CNV_VARIANTS.out.vcf,Channel.value('cnv'))
        ch_versions = ch_versions.mix(CNV_TO_TSV.out.versions)
    }

    if (params.workflow == 'rna'){
        // Annotate gene and transcript tables
        ANNOTATE_EXPRESSION_TABLES(ch_rnaseq_dragen_output, ch_rnaseq_transcript_table)
        ch_versions = ch_versions.mix(ANNOTATE_EXPRESSION_TABLES.out.versions)
    }

    if (params.workflow == 'bsseq'){
        // Need to add a process to convert dragen methylation output to a bed file.
        MAKE_METH_BED (
            DRAGEN_MULTIALIGN.out.dragen_output
        )
        ch_versions = ch_versions.mix(MAKE_METH_BED.out.versions)

        MAKE_METH_BIGWIG (
            DRAGEN_MULTIALIGN.out.dragen_output,
            ch_fasta_reference
        )
        ch_versions = ch_versions.mix(MAKE_METH_BIGWIG.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name    : 'software_versions.yml',
            sort    : true,
            newLine : true
        )
        .set{ ch_collated_versions }
        
    //
    // MODULE: MultiQC
    //
    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_collated_versions
                        .mix(
                            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
                        )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()  // channel: [ path(file) ]
    versions       = ch_versions                  // channel: [ path(file) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
