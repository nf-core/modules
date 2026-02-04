#!/usr/bin/env nextflow

//
// Demultiplex Illumina BCL data using bcl-convert or bcl2fastq
//

include { BCLCONVERT  } from "../../../modules/nf-core/bclconvert/main"
include { BCL2FASTQ   } from "../../../modules/nf-core/bcl2fastq/main"
include { MULTIQC_SAV } from "../../../modules/nf-core/multiqc_sav/main"

workflow BCL_DEMULTIPLEX {
    take:
    ch_flowcell // [[id:"", lane:""], samplesheet.csv, path/to/bcl/files]
    demultiplexer // bclconvert or bcl2fastq

    main:
    ch_fastq = channel.empty()
    ch_reports = channel.empty()
    ch_stats = channel.empty()
    ch_interop = channel.empty()
    ch_logs = channel.empty()
    ch_runinfo = channel.empty()

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
        ch_fastq = ch_fastq.mix(BCLCONVERT.out.fastq)
        ch_interop = ch_interop.mix(BCLCONVERT.out.interop)
        ch_reports = ch_reports.mix(BCLCONVERT.out.reports)
        ch_logs = ch_logs.mix(BCLCONVERT.out.logs)
        ch_runinfo = ch_runinfo.mix(BCLCONVERT.out.runinfo)
    }

    // MODULE: bcl2fastq
    // Demultiplex the bcl files
    if (demultiplexer == "bcl2fastq") {
        BCL2FASTQ(ch_flowcells)
        ch_fastq = ch_fastq.mix(BCL2FASTQ.out.fastq)
        ch_interop = ch_interop.mix(BCL2FASTQ.out.interop)
        ch_reports = ch_reports.mix(BCL2FASTQ.out.reports)
        ch_stats = ch_stats.mix(BCL2FASTQ.out.stats)
        ch_runinfo = ch_runinfo.mix(BCL2FASTQ.out.runinfo)
    }

    // Generate meta for each fastq
    ch_fastq
        .transpose()
        .map { fc_meta, fastq ->
            def meta = [:]
            meta.id = fastq.getSimpleName().toString() - ~/_R[0-9]_001.*$/
            meta.samplename = fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/
            meta.fcid = fc_meta.id
            meta.lane = fc_meta.lane
            // The buffered input stream allows reading directly from cloud storage
            // It will not make a local copy of the file.
            def line = ""
            fastq.withInputStream { fq ->
                def gzipStream = new java.util.zip.GZIPInputStream(fq)
                def decoder = new InputStreamReader(gzipStream, 'ASCII')
                def buffered = new BufferedReader(decoder)
                line = buffered.readLine()
                buffered.close()
            }
            if (line != null && line.startsWith('@')) {
                line = line.substring(1)
                // expected format is like:
                // xx:yy:FLOWCELLID:LANE:... (seven fields)
                def fields = line.split(':')
                // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
                // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
                //def sequencer_serial = fields[0]
                //def run_nubmer       = fields[1]
                def fcid = fields[2]
                def lane = fields[3]
                def index = fields[-1] =~ /[GATC+-]/ ? fields[-1] : ""
                def ID = [fcid, lane].join(".")
                def PU = [fcid, lane, index].findAll().join(".")
                def PL = "ILLUMINA"
                def SM = fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/
                meta.readgroup = [
                    "ID": ID,
                    "SM": SM,
                    "PL": PL,
                    "PU": PU,
                ]
                meta.empty = false
            }
            else {
                println("No reads were found in FASTQ file: ${fastq}")
                meta.readgroup = [:]
                meta.empty = true
            }
            return [meta, fastq]
        }
        .groupTuple(by: [0])
        .map { meta, fastq ->
            meta.single_end = fastq.size() == 1
            return [meta, fastq.flatten()]
        }
        .branch { meta, _fastq ->
            fastq: meta.empty == false
            empty_fastq: meta.empty == true
        }
        .set { ch_fastq_with_meta }

    // MODULE: multiqc_sav
    // Generate MultiQC report for demultiplexing step
    def ch_sav_input = ch_runinfo.join(ch_interop)
    def ch_mqc_input = ch_reports.mix(ch_stats).map { _meta, files -> files.flatten() }
    MULTIQC_SAV(
        ch_sav_input,
        ch_mqc_input,
        [],
        [],
        [],
        [],
        [],
    )

    emit:
    fastq          = ch_fastq_with_meta.fastq
    empty_fastq    = ch_fastq_with_meta.empty_fastq
    reports        = ch_reports
    stats          = ch_stats
    interop        = ch_interop
    logs           = ch_logs
    multiqc_report = MULTIQC_SAV.out.report
    multiqc_data   = MULTIQC_SAV.out.data
    multiqc_plots  = MULTIQC_SAV.out.plots
}
