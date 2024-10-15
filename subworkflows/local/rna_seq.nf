include { DRAGEN_MULTIALIGN as DRAGEN_RNA } from '../../modules/local/dragen_multialign.nf'
include { ANNOTATE_RNASEQ                 } from '../../modules/local/annotate_rnaseq.nf'
include { GET_SIZES_FILE                  } from '../../modules/local/get_sizes_file.nf'
include { GET_STRANDEDNESS                } from '../../modules/local/get_strandedness.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW } from '../../modules/local/bedtools_genomecov.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_RV } from '../../modules/local/bedtools_genomecov.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_US } from '../../modules/local/bedtools_genomecov.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD    } from './bedgraph_bedclip_bedgraphtobigwig.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE    } from './bedgraph_bedclip_bedgraphtobigwig.nf'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_UNSTRANDED } from './bedgraph_bedclip_bedgraphtobigwig.nf'

workflow RNASEQ {
    take:
    input_data
    dragen_inputs

    main:
    ch_versions = Channel.empty()
    ch_dragen_output = Channel.empty()

    DRAGEN_RNA(input_data, dragen_inputs)
    ch_versions = ch_versions.mix(DRAGEN_RNA.out.versions)

    ANNOTATE_RNASEQ(DRAGEN_RNA.out.dragen_output, dragen_inputs)
    ch_versions = ch_versions.mix(ANNOTATE_RNASEQ.out.versions)

    GET_SIZES_FILE(dragen_inputs)
    GET_SIZES_FILE.out.sizes.dump(tag: "sizes")

    GET_STRANDEDNESS(DRAGEN_RNA.out.dragen_output, dragen_inputs)
    bedtools_input =
        GET_STRANDEDNESS.out.dragen_output
        .map{
        meta, dragenoutput, strandedness -> 
        meta.put("strandedness",strandedness)
        return [meta, dragenoutput]
        }
    
    unstranded_input = bedtools_input.filter{ meta, dragen_output -> meta.strandedness == 'unstranded' }
    stranded_input = bedtools_input.filter{ meta, dragen_output -> meta.strandedness != 'unstranded' }

    BEDTOOLS_GENOMECOV_FW(stranded_input, dragen_inputs, 'forward')
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_FW.out.versions)

    BEDTOOLS_GENOMECOV_RV(stranded_input, dragen_inputs, 'reverse')
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_RV.out.versions)

    BEDTOOLS_GENOMECOV_US(unstranded_input, dragen_inputs, 'unstranded')
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_US.out.versions)

    BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD(BEDTOOLS_GENOMECOV_FW.out.bedgraph, GET_SIZES_FILE.out.sizes)

    BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE(BEDTOOLS_GENOMECOV_RV.out.bedgraph, GET_SIZES_FILE.out.sizes)

    BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_UNSTRANDED(BEDTOOLS_GENOMECOV_US.out.bedgraph, GET_SIZES_FILE.out.sizes)
    
    emit:
    dragen_output = DRAGEN_RNA.out.dragen_output
    versions = ch_versions
}