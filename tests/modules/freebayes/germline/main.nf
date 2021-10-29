#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREEBAYES_GERMLINE } from '../../../../modules/freebayes/germline/main.nf' addParams( options: [:] )

workflow test_freebayes {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = []
    samples = []
    populations = []
    cnv = []

    FREEBAYES_GERMLINE ( input, fasta, fai, targets, samples, populations, cnv)
}

workflow test_freebayes_bed {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)]
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai         = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets     = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    samples     = []
    populations = []
    cnv         = []

    FREEBAYES_GERMLINE ( input, fasta, fai, targets, samples, populations, cnv)
}

workflow test_freebayes_cram {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
            ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = []
    samples = []
    populations = []
    cnv = []

    FREEBAYES_GERMLINE ( input, fasta, fai, targets, samples, populations, cnv)
}
