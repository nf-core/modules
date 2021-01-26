#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
def quant_options = [args: '--minAssignedFrags 1']
include { SALMON_INDEX } from '../../../software/salmon/index/main.nf' addParams( options: [:] )
include { SALMON_QUANT } from '../../../software/salmon/quant/main.nf'  addParams( options: quant_options )

workflow test_salmon_index {
    genome_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    transcript_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    SALMON_INDEX ( genome_fasta, transcript_fasta )
}

workflow test_salmon_quant_pe {
    genome_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    transcript_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/covid/GCA_011545545.1_ASM1154554v1_genomic.gtf", checkIfExists: true)
    input = [ [ id:'test', single_end:false ], // meta map
               [ file("${launchDir}/tests/data/fastq/rna/covid/EPI_ISL_486436_1.fastq.gz", checkIfExists: true),
                 file("${launchDir}/tests/data/fastq/rna/covid/EPI_ISL_486436_2.fastq.gz", checkIfExists: true) ] ]

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false)
}

workflow test_salmon_quant_se {
    genome_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_genomic.fna", checkIfExists: true)
    transcript_fasta = file("${launchDir}/tests/data/fasta/covid/GCA_011545545.1_ASM1154554v1_cds_from_genomic.fna", checkIfExists: true)
    gtf = file("${launchDir}/tests/data/gff/covid/GCA_011545545.1_ASM1154554v1_genomic.gtf", checkIfExists: true)
    input = [ [ id:'test', single_end:true ], // meta map
               file("${launchDir}/tests/data/fastq/rna/covid/EPI_ISL_486436_1.fastq.gz", checkIfExists: true) ]

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    SALMON_QUANT( input, SALMON_INDEX.out.index, gtf, transcript_fasta, false)
}
