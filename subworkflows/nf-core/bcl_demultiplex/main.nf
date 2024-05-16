#!/usr/bin/env nextflow

//
// Demultiplex Illumina BCL data using bcl-convert or bcl2fastq
//

include { BCLCONVERT } from "../../../modules/nf-core/bclconvert/main"
include { BCL2FASTQ  } from "../../../modules/nf-core/bcl2fastq/main"

workflow BCL_DEMULTIPLEX {
    take:
        ch_flowcell     // [[id:"", lane:""],samplesheet.csv, path/to/bcl/files]
        demultiplexer   // bclconvert or bcl2fastq

    main:
        ch_versions      = Channel.empty()
        ch_fastq         = Channel.empty()
        ch_reports       = Channel.empty()
        ch_stats         = Channel.empty()
        ch_interop       = Channel.empty()

        // Split flowcells into separate channels containing run as tar and run as path
        // https://nextflow.slack.com/archives/C02T98A23U7/p1650963988498929
        ch_flowcell
            .branch { meta, samplesheet, run ->
                tar: run.toString().endsWith(".tar.gz")
                dir: true
            }.set { ch_flowcells }

        ch_flowcells.tar
            .multiMap { meta, samplesheet, run ->
                samplesheets: [ meta, samplesheet ]
                run_dirs: [ meta, run ]
            }.set { ch_flowcells_tar }

        // Runs when run_dir is a tar archive
        // Re-join the metadata and the untarred run directory with the samplesheet
        ch_flowcells_tar_merged = ch_flowcells_tar
                                    .samplesheets
                                    .join( ch_flowcells_tar.run_dirs )

        // Merge the two channels back together
        ch_flowcells = ch_flowcells.dir.mix(ch_flowcells_tar_merged)

        // MODULE: bclconvert
        // Demultiplex the bcl files
        if (demultiplexer == "bclconvert") {
            BCLCONVERT( ch_flowcells )
            ch_fastq    = ch_fastq.mix(BCLCONVERT.out.fastq)
            ch_interop  = ch_interop.mix(BCLCONVERT.out.interop)
            ch_reports  = ch_reports.mix(BCLCONVERT.out.reports)
            ch_versions = ch_versions.mix(BCLCONVERT.out.versions)
        }

        // MODULE: bcl2fastq
        // Demultiplex the bcl files
        if (demultiplexer == "bcl2fastq") {
            BCL2FASTQ( ch_flowcells )
            ch_fastq    = ch_fastq.mix(BCL2FASTQ.out.fastq)
            ch_interop  = ch_interop.mix(BCL2FASTQ.out.interop)
            ch_reports  = ch_reports.mix(BCL2FASTQ.out.reports)
            ch_stats    = ch_stats.mix(BCL2FASTQ.out.stats)
            ch_versions = ch_versions.mix(BCL2FASTQ.out.versions)
        }

        // Generate meta for each fastq
        ch_fastq_with_meta = ch_fastq
            .map{ fc_meta, fastq -> 
                def meta = [
                    "id": fastq.getSimpleName().toString() - ~/_R[0-9]_001.*$/,
                    "samplename": fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/,
                    "fcid": fc_meta.id,
                    "lane": fc_meta.lane
                ]
                [meta, fastq] 
            }
            .groupTuple(by: [0])
            .map { meta, fastq -> // Add meta.single_end
                meta.single_end = fastq.size() == 1
                return [meta, fastq.flatten()]
            }

    emit:
        fastq    = ch_fastq_with_meta
        reports  = ch_reports
        stats    = ch_stats
        interop  = ch_interop
        versions = ch_versions
}
