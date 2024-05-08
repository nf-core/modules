#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITINTERVALS } from '../../../../../modules/nf-core/gatk4/splitintervals/main.nf'

workflow test_gatk4_splitintervals_bed {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK4_SPLITINTERVALS ( input, fasta, fai, dict)
}

workflow test_gatk4_splitintervals_intervals {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    ]

    GATK4_SPLITINTERVALS ( input, fasta, fai, dict)
}
