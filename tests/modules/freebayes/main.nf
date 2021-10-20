#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREEBAYES } from '../../../modules/freebayes/main.nf' addParams( options: [:] )

workflow test_freebayes {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]
    reference = [file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
                 file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]
    targets = []
    samples = []
    populations = []
    cnv = []

    FREEBAYES ( input, reference, targets, samples, populations, cnv)
}

workflow test_freebayes_bed {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]
    reference = [file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
                 file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]
    targets = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    samples = []
    populations = []
    cnv = []

    FREEBAYES ( input, reference, targets, samples, populations, cnv)
}
