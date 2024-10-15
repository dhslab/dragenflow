include { UCSC_BEDCLIP          } from '../../modules/local/ucsc_bedclip.nf'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/local/ucsc_bedgraphtobigwig.nf'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    //    path: chrom.sizes

    main:
    ch_versions = Channel.empty()

    UCSC_BEDCLIP ( bedgraph, sizes )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.bedgraph, sizes )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

}
