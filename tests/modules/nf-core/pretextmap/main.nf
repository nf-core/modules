#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRETEXTMAP } from '../../../../modules/nf-core/pretextmap/main.nf'

workflow test_pretextmap_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    PRETEXTMAP ( input, [] )
}

workflow test_pretextmap_cram {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    PRETEXTMAP ( input, fasta )
}

workflow test_pretextmap_pairs_gz {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/4dn-dcic/pairix/master/samples/test_4dn.pairs.gz", checkIfExists: true)
    ]

    PRETEXTMAP ( input, [] )
}
