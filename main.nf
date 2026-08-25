#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/dragenflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/dragenflow
    Website: https://nf-co.re/dragenflow
    Slack  : https://nfcore.slack.com/channels/dragenflow
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DRAGENFLOW    } from './workflows/dragenflow.nf'
include { NIRVANA       } from './workflows/nirvana.nf'
include { BUILD_HASH    } from './workflows/build_hash.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_dragenflow_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_dragenflow_pipeline'

//
// WORKFLOW: Run main nf-core/dragenflow analysis pipeline
//
workflow NF_DRAGENFLOW {
    take:
    ch_samplesheet  // channel: [ path(file) ]

    main:
    DRAGENFLOW (ch_samplesheet)

    emit:
    multiqc_report = DRAGENFLOW.out.multiqc_report
    versions       = DRAGENFLOW.out.versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    main:
    ch_versions = Channel.empty()
    ch_multiqc_reports = Channel.value([])

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.demux_outdir,
        params.input
    )
    ch_versions = ch_versions.mix(PIPELINE_INITIALISATION.out.versions)

    if (params.workflow == "download_nirvana") {
        NIRVANA ()

    } else if (params.workflow == "build_hash") {
        BUILD_HASH (params.fasta)

    } else if (params.workflow == "dragenflow") {
        NF_DRAGENFLOW (PIPELINE_INITIALISATION.out.input)
        ch_versions = ch_versions.mix(NF_DRAGENFLOW.out.versions)
        ch_multiqc_reports = NF_DRAGENFLOW.out.multiqc_report
    } else {
        NF_DRAGENFLOW (PIPELINE_INITIALISATION.out.input)
        ch_versions = ch_versions.mix(NF_DRAGENFLOW.out.versions)
        ch_multiqc_reports = NF_DRAGENFLOW.out.multiqc_report
    }

    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ch_multiqc_reports
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
