#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {BAMTOOLS_SPLIT as MARKDUP_P_MERGE_SPLIT            } from "../../modules/bamtools/split/main"
include {BIOBAMBAM_BAMMARKDUPLICATES2 as MARKDUP_P_MARKDUP  } from "../../modules/nf-core/modules/biobambam/bammarkduplicates2/main"
include {SAMTOOLS_MERGE as MARKDUP_P_MERGE                  } from "../../modules/nf-core/modules/samtools/merge/main"

workflow MARKDUP_PARALLEL {
    take:
        ch_input_bam    // [meta, [bam1,bam2,...]]

    main:
        ch_versions = Channel.empty()

        // merge the bams and split per chromosome or interval
        MARKDUP_P_MERGE_SPLIT(ch_input_bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MERGE_SPLIT.out.versions)

        // markduplicates
        MARKDUP_P_MARKDUP(BAMPROCESSING_MERGE_SPLIT.out.bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MARKDUP.out.versions)

        // re-merge bam
        MARKDUP_P_MERGE(BAMPROCESSING_MARKDUP.out.bam)
        ch_versions = ch_versions.mix(BAMPROCESSING_MERGE.out.versions)

        // TODO re-merge metrics
        ch_metrics_merged = Channel.empty()

    emit:
        bam: BAMPROCESSING_MERGE.out.bam    // [meta, bam]
        metrics: ch_metrics_merged          // [meta, metrics]
        versions = ch_versions              // versions
}
