#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREEBAYES } from '../../../modules/freebayes/main.nf'

workflow test_freebayes {
    targets     = []
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              [],
              [],
              targets
            ]
    fasta       = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    samples     = []
    populations = []
    cnv         = []

    FREEBAYES (input, fasta, fai, samples, populations, cnv)
}

workflow test_freebayes_bed {

    targets     = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              [],
              [],
              targets
            ]
    fasta       = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    samples     = []
    populations = []
    cnv         = []

    FREEBAYES (input, fasta, fai, samples, populations, cnv)
}

workflow test_freebayes_cram {

    targets     = []
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
              [],
              [],
              targets
            ]
    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    samples     = []
    populations = []
    cnv         = []

    FREEBAYES (input, fasta, fai, samples, populations, cnv)
}

workflow test_freebayes_somatic {

    targets     = []
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
                targets
            ]
    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    samples     = []
    populations = []
    cnv         = []

    FREEBAYES (input, fasta, fai, samples, populations, cnv)
}

workflow test_freebayes_somatic_cram_intervals {

    targets     = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram_crai'], checkIfExists: true),
                targets
            ]
    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    samples     = []
    populations = []
    cnv         = []

    FREEBAYES (input, fasta, fai, samples, populations, cnv)
}
