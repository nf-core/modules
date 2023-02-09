include { DEFINE_DATASET_FROM_ANNO  } from '../../../modules/nf-core/aberrantsplicing/precomp/definedatasetfromanno'

workflow ABERRANT_SPLICING {
    take:
        params

    main:
        ch_versions = Channel.empty()

        Channel.from(params.gtf.split(", "))
            .map{iM -> iM.split(":")}
            .set{addProc}
        groups = params.aberrantslicing.groups.split(", ")
        sampleAnn = params.sample_annotation_relative

        // Counting part
        Channel.fromPath(params.bams).set{ch_bam}
        DEFINE_DATASET_FROM_ANNO(sampleAnn, ch_bam.collect())
}
