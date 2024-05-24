include { SAMPLESHEET_CHECK              } from '../../modules/local/samplesheet_check.nf'
include { GATHER_FASTQS                  } from '../../subworkflows/local/gather_fastqs.nf'

def create_master_samplesheet(LinkedHashMap row) {

    // create meta map
    def meta = [:]

    def meta = [:]
    meta.id             = row.id
    meta.uid            = row.uid ?: null
    meta.sample_type    = row.sample_type ?: null
    meta.sample_id      = row.sample_id ?: null
    meta.assay          = row.assay ?: null

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
    count = 0
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


workflow SOMATIC_INPUT_CHECK {
    take:
    master_samplesheet
    data_path

    main:

    ch_mastersheet        = Channel.empty()
    ch_input_data         = Channel.empty()
    ch_dragen_outputs     = Channel.empty()

    // Runs a python script that parses the sample sheet and adds key metadata, 
    // including index sequences, flowcell, and lane. If fastq_list.csv files are passed,
    // these are parsed also and read1/read2 pairs are returned. 
    // Also, if tumor/normal and case ids are supplied, mark the 
    // rows in the sheet as tumor or normal and assign a case id.
    // The output is then channelified.
    SAMPLESHEET_CHECK ( master_samplesheet, data_path )
    .csv
    .splitCsv ( header:true, sep:',' )
    .map { create_master_samplesheet(it) }
    .map { meta, inputs -> 
        if (meta.id == params.tumorid){
            meta.sample_type = 'tumor'
            meta.uid = params.id
        } else if (meta.id == params.normalid){
            meta.sample_type = 'normal'
            meta.uid = params.id
        }   
        return [ meta, inputs ]
    }
    .filter { it[0].sample_type != null || it[0].dragen_path != null }
    .groupTuple()
    .map { meta, inputs -> 
        def new_inputs = merge_maps(inputs)
        meta.count = new_inputs.count
        [ meta, new_inputs ] 
    }
    .set { ch_mastersheet }

    // we need fastq lists for somatic mode. Cant align 2 cram/bam files. 
    // Gather fastqs returns [ meta, fastq_list] 
    GATHER_FASTQS ( ch_mastersheet )

    ch_fastq_lists
    .map { meta, fqlist -> 
        def key = groupKey(meta, meta.count)
        [ key, fqlist ]
    }
    .groupTuple() | CONCATENATE_FASTQLISTS // this is done with a process because otherwise it could wait for the group to fill */

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
        [ meta, files ]
    }
    .set { ch_align_inputs }
    
    // get read info for fastq_list file generation
    ch_mastersheet
    .map { meta, files -> 
        def new_meta = [:]
        new_meta['id'] = meta.uid
        new_meta['assay'] = meta.assay
        tumor_id = ""
        normal_id = ""
        if (meta.sample_type == 'tumor'){
            tumor_id = meta.id
        } else if (meta.sample_type == 'normal'){
            normal_id = meta.id
        }

        [ new_meta, tumor_id, normal_id, files ]
    }
    .groupTuple(by:0, size=2)
    .map { meta, tumor, normal, files -> 
            new_meta = meta.subMap('id', 'assay')
            new_meta['tumor'] = tumor.findAll { it != '' }.unique()[0]
            new_meta['normal'] = normal.findAll { it != '' }.unique()[0]

            return [ new_meta, files.flatten() ]
    }
    .filter { it[0].tumor != "" && it[0].normal != "" }
    .set { ch_fastqs }

    ch_mastersheet
    .map { meta, inputs -> 
        if (inputs.dragen_path != null){
            def new_meta = meta.subMap('id','assay')
            return [ new_meta, file(inputs.dragen_path).listFiles() ]
        }
    }
    .set { ch_dragen_outputs }

    emit:
    dragen_outputs = ch_dragen_outputs
    input_data = ch_input_data
}
