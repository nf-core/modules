#!/usr/bin/env nextflow

//
// Demultiplex Illumina BCL data using bcl-convert or bcl2fastq
//

include { BCLCONVERT } from "../../../modules/nf-core/bclconvert/main"
include { BCL2FASTQ  } from "../../../modules/nf-core/bcl2fastq/main"

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

def generateReadgroupBCLCONVERT(ch_fastq_list_csv, ch_fastq) {
    return ch_fastq_list_csv
        .join(ch_fastq, by: [0])
        .map { meta, csv_file, fastq_list ->
            def meta_fastq = []
            csv_file
                .splitCsv(header: true)
                .each { row ->
                    // Create the readgroup tuple
                    // RGID,RGSM,RGLB,Lane,Read1File,Read2File
                    def rg = [:]
                    // row.RGID is index1.index2.lane
                    rg.ID = row.RGID
                    // RGPU is a custom column in the samplesheet containing the flowcell ID
                    rg.PU = row.RGPU ? row.RGPU : meta.id + "." + row.Lane
                    rg.SM = row.RGSM
                    rg.LB = row.RGLB ? row.RGLB : ""
                    rg.PL = "ILLUMINA"

                    // dereference the fastq files in the csv
                    def fastq1 = fastq_list.find { fq -> file(fq).name == file(row.Read1File).name }
                    def fastq2 = row.Read2File ? fastq_list.find { fq -> file(fq).name == file(row.Read2File).name } : null

                    // set fastq metadata
                    def new_meta = meta + [id: fastq1.getSimpleName().toString() - ~/_R[0-9]_001.*$/, readgroup: rg, single_end: !fastq2]

                    meta_fastq << [new_meta, fastq2 ? [fastq1, fastq2] : [fastq1]]
                }
            return meta_fastq
        }
        .flatMap()
}

def generateReadgroupBCL2FASTQ(ch_fastq) {
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
                def ID = [index, lane].join(".")
                def LB = ""
                def PL = "ILLUMINA"
                def PU = [fcid, lane].findAll().join(".")
                def SM = fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/
                meta.readgroup = ["ID": ID, "SM": SM, "PL": PL, "PU": PU, "LB": LB]
            }
            else {
                println("No reads were found in FASTQ file: ${fastq}")
                meta.readgroup = [:]
            }
            return [meta, fastq]
        }
        .groupTuple(by: [0])
        .map { meta, fastq ->
            meta.single_end = fastq.size() == 1
            return [meta, fastq.flatten()]
        }
}
