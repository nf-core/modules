#!/usr/bin/env nextflow

//
// Demultiplex Illumina BCL data using bcl-convert or bcl2fastq
//

include { BCLCONVERT                                       } from "../../../modules/nf-core/bclconvert/main"
include { generateReadgroup as generateReadgroupBCLCONVERT } from "../../../modules/nf-core/bclconvert/main"
include { BCL2FASTQ                                        } from "../../../modules/nf-core/bcl2fastq/main"
include { generateReadgroup as generateReadgroupBCL2FASTQ  } from "../../../modules/nf-core/bcl2fastq/main"

workflow BCL_DEMULTIPLEX {
    take:
    ch_flowcell // [[id:"", lane:""], samplesheet.csv, path/to/bcl/files]
    demultiplexer // bclconvert or bcl2fastq

    main:
    ch_fastq_with_meta = channel.empty()
    ch_reports = channel.empty()
    ch_stats = channel.empty()
    ch_interop = channel.empty()
    ch_logs = channel.empty()

    // Split flowcells into separate channels containing run as tar and run as path
    // https://nextflow.slack.com/archives/C02T98A23U7/p1650963988498929
    ch_flowcell
        .branch { _meta, _samplesheet, run ->
            tar: run.toString().endsWith(".tar.gz")
            dir: true
        }
        .set { ch_flowcells }

    ch_flowcells.tar
        .multiMap { meta, samplesheet, run ->
            samplesheets: [meta, samplesheet]
            run_dirs: [meta, run]
        }
        .set { ch_flowcells_tar }

    // Runs when run_dir is a tar archive
    // Re-join the metadata and the untarred run directory with the samplesheet
    ch_flowcells_tar_merged = ch_flowcells_tar.samplesheets.join(ch_flowcells_tar.run_dirs)

    // Merge the two channels back together
    ch_flowcells = ch_flowcells.dir.mix(ch_flowcells_tar_merged)

    // MODULE: bclconvert
    // Demultiplex the bcl files
    if (demultiplexer == "bclconvert") {
        BCLCONVERT(ch_flowcells)
        ch_interop = ch_interop.mix(BCLCONVERT.out.interop)
        ch_reports = ch_reports.mix(BCLCONVERT.out.reports)
        ch_logs = ch_logs.mix(BCLCONVERT.out.logs)
        ch_fastq_with_meta = ch_fastq_with_meta.mix(
            generateReadgroupBCLCONVERT(
                BCLCONVERT.out.reports.map { meta, reports ->
                    return [meta, reports.find { report -> report.name == "fastq_list.csv" }]
                },
                BCLCONVERT.out.fastq,
            )
        )
    }

    // MODULE: bcl2fastq
    // Demultiplex the bcl files
    if (demultiplexer == "bcl2fastq") {
        BCL2FASTQ(ch_flowcells)
        ch_interop = ch_interop.mix(BCL2FASTQ.out.interop)
        ch_reports = ch_reports.mix(BCL2FASTQ.out.reports)
        ch_stats = ch_stats.mix(BCL2FASTQ.out.stats)

        ch_fastq_with_meta = ch_fastq_with_meta.mix(
            generateReadgroupBCL2FASTQ(
                BCL2FASTQ.out.fastq
            )
        )
    }

    // extract empty fastq files from channel
    ch_fastq = ch_fastq_with_meta.branch { meta, fastq ->
        empty: fastq.any { fq -> file(fq).size() < 30 }
            return [meta, fastq]
        fastq: true
            return [meta, fastq]
    }

    emit:
    fastq       = ch_fastq.fastq
    empty_fastq = ch_fastq.empty
    reports     = ch_reports
    stats       = ch_stats
    interop     = ch_interop
    logs        = ch_logs
}
