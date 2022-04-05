#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCLCONVERT } from '../../../modules/bclconvert/main.nf'

process STUB_BCLCONVERT {
    output:
    path "*.fastq.gz"               ,emit: fastq
    path "Reports/*.{csv,xml,bin}"  ,emit: reports
    path "Logs/*.{log,txt}"         ,emit: logs
    path "versions.yml"             ,emit: versions

    stub:
    """
    touch sample1_S1_L001_R1_001.fastq.gz
    touch sample1_S1_L001_R2_001.fastq.gz
    touch sample1_S1_L002_R1_001.fastq.gz
    touch sample1_S1_L002_R2_001.fastq.gz
    touch sample2_S2_L001_R1_001.fastq.gz
    touch sample2_S2_L001_R2_001.fastq.gz
    touch sample2_S2_L002_R1_001.fastq.gz
    touch sample2_S2_L002_R2_001.fastq.gz

    mkdir Reports
    touch Reports/Adapter_Metrics.csv
    touch Reports/Demultiplex_Stats.csv
    touch Reports/fastq_list.csv
    touch Reports/Index_Hopping_Counts.csv
    touch Reports/IndexMetricsOut.bin
    touch Reports/Quality_Metrics.csv
    touch Reports/RunInfo.xml
    touch Reports/SampleSheet.csv
    touch Reports/Top_Unknown_Barcodes.csv

    mkdir Logs
    touch Logs/Errors.log
    touch Logs/FastqComplete.txt
    touch Logs/Info.log
    touch Logs/Warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: 00.000.000.3.9.3
    END_VERSIONS
    """
}

workflow test_bclconvert {
    STUB_BCLCONVERT ()
}
