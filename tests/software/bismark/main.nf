#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISMARK_GENOME_PREPARATION        } from '../../../software/bismark/genome_preparation/main.nf' addParams( options: [:] )
include { BISMARK_ALIGN as BISMARK_ALIGN_SE } from '../../../software/bismark/align/main.nf'  addParams( options: [ publish_dir:'test_single_end' ] )
include { BISMARK_ALIGN as BISMARK_ALIGN_PE } from '../../../software/bismark/align/main.nf'  addParams( options: [ publish_dir:'test_paired_end' ] )
include { BISMARK_DEDUPLICATE               } from '../../../software/bismark/deduplicate/main.nf'  addParams( options: [:] )
include { BISMARK_EXTRACT                   } from '../../../software/bismark/extract/main.nf'  addParams( options: [:] )
include { BISMARK_REPORT                    } from '../../../software/bismark/report/main.nf'  addParams( options: [:] )
include { BISMARK_SUMMARY                   } from '../../../software/bismark/summary/main.nf'  addParams( options: [:] )

workflow test_bismark_genome_preparation {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) ]

    BISMARK_GENOME_PREPARATION ( input )
}

/*
 * Test with single-end data
 */
workflow test_bismark_align_single_end {

    def input = []
    def index = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true) ] ]
    index = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/index/E_coli/bismark/BismarkIndex", checkIfExists: true) ] ]

    BISMARK_ALIGN_SE (
        input,
        index
    )
}

/*
 * Test with paired-end data
 */
workflow test_bismark_align_paired_end {

    def input = []
    def index = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/methylated_dna/Ecoli_10K_methylated_R2.fastq.gz", checkIfExists: true) ] ]
    index = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/index/E_coli/bismark/BismarkIndex", checkIfExists: true) ] ]

    BISMARK_ALIGN_PE (
        input,
        index
    )
}

workflow test_bismark_deduplicate {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.bam", checkIfExists: true) ] ]

    BISMARK_DEDUPLICATE ( input )
}

workflow test_bismark_extract {

    def input = []
    def index = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.bam", checkIfExists: true) ] ]
    index = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/index/E_coli/bismark/BismarkIndex", checkIfExists: true) ] ]

    BISMARK_EXTRACT (
        input,
        index
    )
}

workflow test_bismark_report {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_PE_report.txt", checkIfExists: true),
              file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.deduplication_report.txt", checkIfExists: true),
              file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe_splitting_report.txt", checkIfExists: true),
              file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.M-bias.txt", checkIfExists: true) ]

    BISMARK_REPORT ( input )
}

workflow test_bismark_summary {

    def bam = file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.bam", checkIfExists: true)
    def align = file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_PE_report.txt", checkIfExists: true)
    def dedup = file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.deduplication_report.txt", checkIfExists: true)
    def split = file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe_splitting_report.txt", checkIfExists: true)
    def mbias = file("${launchDir}/tests/data/bismark/Ecoli_10K_methylated_R1_bismark_bt2_pe.M-bias.txt", checkIfExists: true)

    BISMARK_SUMMARY (
        bam,
        align,
        dedup,
        split,
        mbias
    )
}
