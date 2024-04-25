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

WorkflowDragenmultiworkflow.initialise(params, log)

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
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK          } from '../subworkflows/local/input_check.nf'
include { SOMATIC_INPUT_CHECK  } from '../subworkflows/local/somatic_input_check.nf'
include { TUMOR_NORMAL         } from '../subworkflows/local/tumor_normal.nf'
include { TUMOR                } from '../subworkflows/local/tumor.nf'
//include { GERMLINE             } from '../subworkflows/local/germline.nf'
include { RNASEQ               } from '../subworkflows/local/rna_seq.nf'
include { METHYLATION          } from '../subworkflows/local/methylation.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This function 'stages' a set of files defined by a map of key:filepath pairs.
// It returns a tuple: a map of key:filename pairs and list of file paths.
// This can be used to generate a value Channel that can be used as input to a process
// that accepts a tuple val(map), path("*") so map.key refers to the appropriate linked file.
def stageFileset(Map filePathMap) {
    def basePathMap = [:]
    def filePathsList = []

    filePathMap.each { key, value ->
        if (value != null) {
            def filepath = file(value)
            if (filepath.exists()) {
                // Add basename and key to the map
                basePathMap[key] = value.split('/')[-1]
                // Add file path to the list
                filePathsList << filepath
            } else {
                println "Warning: File at '${value}' for key '${key}' does not exist."
            }
        }
    }
    return [basePathMap, filePathsList]
}

// If MGI samplesheet is used, we need to set the 
// data path because only files are given. This sets the 
// data path to the samplesheet directory, or the data_path parameter.
def data_path = ""
def mastersheet = params.mastersheet
if (params.mgi == true) {
    data_path = new File(params.mastersheet).parentFile.absolutePath
} else if (params.data_path != null){
    data_path  = params.data_path
}

// Info required for completion email and summary
def multiqc_report = []

workflow DRAGENMULTIWORKFLOW {

    ch_versions     = Channel.empty()
    ch_input_data   = Channel.empty()

    if (params.workflow == 'somatic') {

        SOMATIC_INPUT_CHECK(Channel.fromPath(mastersheet), data_path)
        ch_input_data = ch_input_data.mix(SOMATIC_INPUT_CHECK.out.input_data)
        ch_versions = ch_versions.mix(SOMATIC_INPUT_CHECK.out.versions)

        TUMOR_NORMAL(ch_input_data, ch_dragen_inputs)
        ch_versions = ch_versions.mix(TUMOR_NORMAL.out.ch_versions)
    
    } else {

        INPUT_CHECK(Channel.fromPath(mastersheet), data_path)
        ch_input_data = ch_input_data.mix(INPUT_CHECK.out.input_data)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)        

        if (params.workflow == '5mc') {

            // Stage Dragen input files
            params.dragen_inputs.reference = params.dragen_inputs.methylation_reference
            params.dragen_inputs.methylation_reference = null
            ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))

            METHYLATION(ch_input_data, ch_dragen_inputs)
            ch_versions = ch_versions.mix(METHYLATION.out.ch_versions)

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
                ch_versions = ch_versions.mix(RNASEQ.out.ch_versions)
            }

            if (params.workflow == 'germline') {
                //GERMLINE(ch_input_data, ch_dragen_inputs)
                //ch_versions = ch_versions.mix(GERMLINE.out.ch_versions)
            }

            if (params.workflow == 'tumor') {
                TUMOR(ch_input_data, ch_dragen_inputs)
                ch_versions = ch_versions.mix(TUMOR.out.ch_versions)
            }
        }
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDragenmultiworkflow.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDragenmultiworkflow.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
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
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
