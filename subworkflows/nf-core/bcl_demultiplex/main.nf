#!/usr/bin/env nextflow

//
// Demultiplex Illumina BCL data using bcl-convert or bcl2fastq
//

include { BCLCONVERT } from "../../../modules/nf-core/bclconvert/main"
include { BCL2FASTQ  } from "../../../modules/nf-core/bcl2fastq/main"

// Define the log file path before the workflow starts
def logFile = file("${params.outdir}/invalid_fastqs.log")

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
        ch_fastq_with_meta = generate_fastq_meta(ch_fastq, logFile)

    emit:
        fastq    = ch_fastq_with_meta
        reports  = ch_reports
        stats    = ch_stats
        interop  = ch_interop
        versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This function appends a given text to a specified log file.
// If the log file does not exist, it creates a new one.
def appendToLogFile(String text, Path logFile) {
    // Convert the text to String if it's a GString
    String textToWrite = text.toString()
    // Append text to the log file
    logFile.withWriterAppend { writer ->
        writer.println(textToWrite) // Appends the text to the file with a new line
    }
}

// Add meta values to fastq channel and skip invalid FASTQ files
def generate_fastq_meta(ch_reads, logFile) {
    // Create a tuple with the meta.id and the fastq
    ch_reads.transpose().map { fc_meta, fastq ->
        // Check if the FASTQ file is empty or has invalid content
        def isValid = fastq.withInputStream { is ->
            new java.util.zip.GZIPInputStream(is).withReader('ASCII') { reader ->
                def line = reader.readLine()
                line != null && line.startsWith('@')
            }
        }

        def meta = null
        if (isValid) {
            meta = [
                "id": fastq.getSimpleName().toString() - ~/_R[0-9]_001.*$/,
                "samplename": fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/,
                "readgroup": [:],
                "fcid": fc_meta.id,
                "lane": fc_meta.lane
            ]
            meta.readgroup = readgroup_from_fastq(fastq)
            meta.readgroup.SM = meta.samplename
        } else {
                appendToLogFile(
                    "Empty or invalid FASTQ file: ${fastq}",
                    logFile
                )
            fastq = null
        }

        return [meta, fastq]
    }.filter { it[0] != null }
    // Group by meta.id for PE samples
    .groupTuple(by: [0])
    // Add meta.single_end
    .map { meta, fastq ->
        if (meta != null) {
            meta.single_end = fastq.size() == 1
        }
        return [meta, fastq.flatten()]
    }
}

// https://github.com/nf-core/sarek/blob/7ba61bde8e4f3b1932118993c766ed33b5da465e/workflows/sarek.nf#L1014-L1040
def readgroup_from_fastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)

    def line

    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    def rg = [:]

    // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
    sequencer_serial = fields[0]
    run_nubmer       = fields[1]
    fcid             = fields[2]
    lane             = fields[3]
    index            = fields[-1] =~ /[GATC+-]/ ? fields[-1] : ""

    rg.ID = [fcid,lane].join(".")
    rg.PU = [fcid, lane, index].findAll().join(".")
    rg.PL = "ILLUMINA"

    return rg
}
