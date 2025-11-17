/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDragenflow.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS
//
include { GATHER_ALIGNMENT_SAMPLES  } from '../subworkflows/local/gather_alignment_samples.nf'
include { CAT_FASTQLISTS            } from '../subworkflows/local/cat_fastqlists.nf'
include { ALIGN                     } from '../subworkflows/local/align.nf'
include { GERMLINE                  } from '../subworkflows/local/germline.nf'
include { METHYLATION               } from '../subworkflows/local/methylation.nf'
include { RNASEQ                    } from '../subworkflows/local/rna_seq.nf'
include { SOMATIC                   } from '../subworkflows/local/somatic.nf'
include { TUMOR                     } from '../subworkflows/local/tumor.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { PARSE_INPUT_SAMPLESHEET     } from '../modules/local/parse_input_samplesheet'
include { SAMPLESHEET_CHECK           } from '../modules/local/samplesheet_check.nf'
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
ch_dragen_hotspots = params.dragen_hotspots
    ? Channel.fromPath("${params.dragen_hotspots}*", checkIfExists: true).collect()
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

    if (params.workflow == 'somatic') {

        // for somatic workflow, need to assemble tumor and normal data
        ch_alignment_samples
        .branch {
            tumor: it[0].sample_type == 'tumor'
            normal: it[0].sample_type == 'normal'
        }.set { ch_samples }

        ch_somatic_input = Channel.empty()

        ch_somatic_input = ch_samples.tumor // first get tumor samples
        .map { meta, reads, fastqlist, alignment_files ->
            [ meta.individual_id, meta, reads, fastqlist ]
        }
        .cross( // and join with normal samples
            ch_samples.normal
            .map { meta, reads, fastqlist, alignment_files ->
                [ meta.individual_id, meta, reads, fastqlist ]
            }
        ).map { tumor, normal -> [ tumor[1], tumor[2], tumor[3], normal[1], normal[2], normal[3] ] }
        // on joined samples, set tumor and normal id to meta and combine reads.
        .map { tumor_meta, tumor_reads, tumor_fastqlist, normal_meta, normal_reads, normal_fastqlist ->
            def new_meta = [:]
            new_meta['id'] = tumor_meta['id']
            new_meta['sex'] = tumor_meta['sex']
            new_meta['tumor_id'] = tumor_meta['id']
            new_meta['normal_id'] = normal_meta['id']
            reads = tumor_reads + normal_reads
            [ new_meta, reads, [ tumor_fastqlist, normal_fastqlist ] ]
        }
        // now separate reads and fastqlists
        .multiMap { meta, reads, fastqlists -> 
            reads: [ meta, reads ]
            fastqlist: [ meta.id, meta, fastqlists ]
        }
        
        ch_alignment_input = ch_somatic_input.reads
            .join(
                CAT_FASTQLISTS (
                    ch_somatic_input.fastqlist
                )
                .fastqlist
            ).dump(tag: 'somatic_alignment_input', pretty: true)
        /*ch_inputs = ch_somatic_input.reads
            .join(
                ch_somatic_input.fastqlist
                .map { meta, tumor_fastqlist, normal_fastqlist -> [ meta.id, meta ] }
                .join(
                    ch_somatic_input.fastqlist
                    .map { meta, tumor_fastqlist, normal_fastqlist ->
                            def tumordata = parseFastqList(tumor_fastqlist)
                            def normaldata = parseFastqList(normal_fastqlist)
                            if (tumordata && normaldata) {
                                def header = tumordata[0].keySet().join(',')
                                def tumorcontent = tumordata.collect { it.values().join(',') }.join('\n')
                                def normalcontent = normaldata.collect { it.values().join(',') }.join('\n')
                                def content = tumorcontent + '\n' + normalcontent                            
                                [ meta, header, content ]
                            }
                    }
                    .collectFile { meta, header, content -> [ "${meta.id}.fastq_list.csv", header + '\n' + content ], keepHeader: true }
                    .map { filename -> [ filename.toString.split('/')[-1].split('.')[0], filename ] }
                )
                .map { id, meta, fastqfile -> [ meta, fastqfile ] }
            )
                
        ch_inputs.dump(pretty: true)
    */
    }
        /*
        SOMATIC(ch_somatic_input, ch_dragen_inputs)
        ch_versions = ch_versions.mix(SOMATIC.out.versions)

    } else {

        if (params.workflow == '5mc') {

            // Stage Dragen input files
            params.dragen_inputs.reference = params.dragen_inputs.methylation_reference
            params.dragen_inputs.methylation_reference = null
            ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))

            METHYLATION(ch_input_data, ch_dragen_inputs)
            ch_versions = ch_versions.mix(METHYLATION.out.versions)

        } else {

            params.dragen_inputs.methylation_reference = null
            if (params.target_bed_file != null){
                params.dragen_inputs.target_bed_file = params.target_bed_file
            }
            if (params.hotspot_vcf != null){
                params.dragen_inputs.hotspot_vcf = params.hotspot_vcf
                params.dragen_inputs.hotspot_vcf_index = params.hotspot_vcf_index
            }

            ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))

            if (params.workflow == 'rna') {
                RNASEQ(ch_input_data, ch_dragen_inputs)
                ch_versions = ch_versions.mix(RNASEQ.out.versions)
            }

            if (params.workflow == 'tumor') {
                ch_input_data.view()
                TUMOR(ch_input_data, ch_dragen_inputs)
                ch_versions = ch_versions.mix(TUMOR.out.versions)
            }

            if (params.workflow == 'align' || (params.workflow == 'idtumis' && params.target_bed_file != null)) {
                ALIGN(ch_input_data, ch_dragen_inputs)
                ch_versions = ch_versions.mix(ALIGN.out.versions)
            }

            if (params.workflow == 'germline'){ 
                GERMLINE(ch_input_data, ch_dragen_inputs)
                ch_versions = ch_versions.mix(GERMLINE.out.versions)
            }

        }
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDragenflow.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDragenflow.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    */

    emit:
    multiqc_report = Channel.empty() //MULTIQC.out.report.toList()  // channel: [ path(file) ]
    versions       = ch_versions                  // channel: [ path(file) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
