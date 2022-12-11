#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTWGSMETRICS } from '../../../../../modules/nf-core/picard/collectwgsmetrics/main.nf'

//Test without an interval list file
workflow test_picard_collectwgsmetrics {
    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            ]
    fasta = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    PICARD_COLLECTWGSMETRICS ( input, fasta, fai,  [] )
}

//Test with an interval list file
workflow test_picard_collectwgsmetrics_with_interval {
    input        =  [ [ id:'test', single_end:false ], // meta map
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                    []
                    ]
    fasta = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    intervallist = file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true)

    PICARD_COLLECTWGSMETRICS ( input, fasta, fai, intervallist )
}
