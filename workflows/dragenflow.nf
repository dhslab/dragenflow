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
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SAMPLESHEET_CHECK         } from '../modules/local/samplesheet_check.nf'
include { GATHER_FASTQS             } from '../subworkflows/local/gather_fastqs.nf'
include { CONCATENATE_FASTQLISTS    } from '../modules/local/concatenate_fastqlists.nf'
include { SOMATIC                   } from '../subworkflows/local/somatic.nf'
include { TUMOR                     } from '../subworkflows/local/tumor.nf'
include { ALIGN                     } from '../subworkflows/local/align.nf'
include { RNASEQ                    } from '../subworkflows/local/rna_seq.nf'
include { METHYLATION               } from '../subworkflows/local/methylation.nf'
//include { GERMLINE             } from '../subworkflows/local/germline.nf'

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

def create_samplesheet(LinkedHashMap row) {
    // create meta map
    def meta = [:]

    meta.id             = row.id ?: null
    meta.assay          = row.assay ?: null
    meta.uid            = row.uid ?: null  
    meta.sample_type    = row.sample_type ?: null
    meta.sample_id      = row.sample_id ?: null
    meta.sex            = row.sex ?: null

    def obj = [:]

    obj.indexes = row.lane && row.i7index && row.i5index ? [ lane:row.lane, i7index:row.i7index, i5index:row.i5index ] : null
    obj.fastq_list = row.fastq_list ? row.fastq_list : null
    obj.demux_path = row.demux_path ? row.demux_path : null
    obj.reads = row.read1 && row.read2 ? [ read1:row.read1, read2:row.read2 ] : null
    obj.cram = row.cram ? row.cram : null
    obj.dragen_path = row.dragen_path ? row.dragen_path : null

    return [ meta, obj ] //indexes, fastq_list, demux_path, reads, cram, dragen_path ]
}

// Function to merge multiple LinkedHashMaps into lists by key using the '<<' operator
LinkedHashMap merge_maps(List<LinkedHashMap> maps) {
    LinkedHashMap result = new LinkedHashMap()
    def count = 0
    maps.each { map ->
        map.each { key, value ->
            if (value != null){        
                if (result.containsKey(key) && result[key] != null) {
                    result[key] << value
                } else {
                    result[key] = [value]
                    count++
                }
            }
        }
    }
    result['count'] = count
    return result
}

def parseCSV(filePath) {
    List<Map<String, String>> data = []

    // Read the file and split each line
    filePath.withReader { reader ->
        // Read the header line to use as keys for the map
        def headers = reader.readLine().split(',')

        // Process each subsequent line
        reader.splitEachLine(',') { values ->
            def row = [:]
            headers.eachWithIndex { header, i ->
                row[header.trim()] = values[i].trim()
            }
            data.add(row)
        }
    }

    return data
}

// If MGI samplesheet is used, we need to set the 
// data path because only files are given. This sets the 
// data path to the samplesheet directory, or the data_path parameter.
def data_path = ""
if (params.mgi == true) {
    data_path = new File(params.input).parentFile.absolutePath
} else if (params.data_path != null){
    data_path  = params.data_path
}

// Info required for completion email and summary
def multiqc_report = []

workflow DRAGENFLOW {

    ch_versions     = Channel.empty()
    ch_mastersheet  = Channel.empty()
    ch_input_data   = Channel.empty()

    // Parse mastersheet and get meta hash and inputs.
    SAMPLESHEET_CHECK(Channel.fromPath(params.input), data_path)
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    SAMPLESHEET_CHECK.out.csv
    .splitCsv ( header:true, sep:',', quote:'"' )
    .map { create_samplesheet(it) }
    .groupTuple()
    .map { meta, inputs -> 
        def new_inputs = merge_maps(inputs)
        meta.count = new_inputs.count
        [ meta, new_inputs ] 
    }
    .branch {
        crams: it[0].count == 1 && it[1].cram != null && params.workflow != 'somatic'
        fastqs: it[1].cram == null || params.workflow == 'somatic'
        other: true
    }.set { ch_mastersheet }

    // get crams. first branch on whether there are multiple crams or a single one.
    // for now do not handle multiple crams.
    ch_mastersheet.crams
    .branch { 
        single: it[1].cram.size() == 1
        multiple: it[1].cram.size() > 1
    }.set { ch_aligned_crams }    

    ch_aligned_crams.single
    .map { meta, inputs -> 
        def new_meta = meta.subMap('id', 'uid','sex','sample_type','sample_id','assay')
        if (inputs.cram[0].split('\\.')[-1] == 'cram') {
            new_meta.cram = file(inputs.cram[0]).getName()
            return [ new_meta, 'cram', [ file(inputs.cram[0]), file(inputs.cram[0] + '.crai', checkIfExists: true)  ] ]
        } else if (inputs.cram[0].split('\\.')[-1] == 'bam') {
            new_meta.bam = file(inputs.cram[0]).getName()
            return [ new_meta, 'bam', [ file(inputs.cram[0]), file(inputs.cram[0] + '.bai', checkIfExists: true)  ] ]
        }
    }
    .mix(ch_input_data)
    .set { ch_input_data }
    
    ch_mastersheet.fastqs.dump(pretty: true)

    //
    // get fastqs
    //
    GATHER_FASTQS(ch_mastersheet.fastqs)
    ch_versions = ch_versions.mix(GATHER_FASTQS.out.versions)

    GATHER_FASTQS.out.fastqs
    .map { meta, fqlist -> 
        def key = groupKey(meta, meta.count)
        [ key, fqlist ]
    }
    .groupTuple() | CONCATENATE_FASTQLISTS // this is done with a process because otherwise it could wait for the group to fill
   
    // this takes the fastq lists and adds the read1/read2 files
    CONCATENATE_FASTQLISTS.out.fastqlist
    .map { meta, fastqlist ->
        def files = [ file(fastqlist) ]
        def listdata = parseCSV(file(fastqlist))
        listdata.each { row ->
            if (row.containsKey('Read1File')) {
                files.add(file(row['Read1File']))
            }
            if (row.containsKey('Read2File')) {
                files.add(file(row['Read2File']))
            }
        }
        [ meta, 'fastq', files ]
    }
    .mix(ch_input_data)
    .set { ch_input_data }

    ch_input_data.dump(pretty: true)

    if (params.workflow == 'somatic') {

        // for somatic workflow, need to assemble tumor and normal data
        // note that for the somatic workflow, the input type is only fastq
        ch_somatic_input = Channel.empty()

        ch_input_data
        .map { meta, type, files -> 
            def new_meta = [:]
            new_meta['id'] = meta.uid
            new_meta['assay'] = meta.assay
            def tumor_id = ""
            def normal_id = ""
            if (meta.sample_type == 'tumor'){
                tumor_id = meta.id
            } else if (meta.sample_type == 'normal'){
                normal_id = meta.id
            }
            [ new_meta, tumor_id, normal_id, files ]
        }
        .groupTuple(by:0)
        .filter { it[1].findAll { it != '' }.unique().size() == 1 && it[2].findAll { it != '' }.unique().size() == 1 }
        .map { meta, tumor, normal, files -> 
                new_meta = meta.subMap('id', 'assay')
                new_meta['tumor'] = tumor.findAll { it != '' }.unique()[0]
                new_meta['normal'] = normal.findAll { it != '' }.unique()[0]

                return [ new_meta, "fastq", files.flatten() ]
        }
        .filter { it[0].tumor != "" && it[0].normal != "" }
        .set { ch_somatic_input }
            
        params.dragen_inputs.methylation_reference = null
        if (params.target_bed_file != null){
            params.dragen_inputs.target_bed_file = params.target_bed_file
        }
        if (params.hotspot_vcf != null){
            params.dragen_inputs.hotspot_vcf = params.hotspot_vcf
            params.dragen_inputs.hotspot_vcf_index = params.hotspot_vcf_index
        }

        ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))

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

            if (params.workflow == 'germline') {
                //GERMLINE(ch_input_data, ch_dragen_inputs)
                //ch_versions = ch_versions.mix(GERMLINE.out.ch_versions)
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
