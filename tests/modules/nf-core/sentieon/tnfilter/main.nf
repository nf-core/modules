#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_TNFILTER } from '../../../../../modules/nf-core/sentieon/tnfilter/main.nf'

workflow test_sentieon_tnfilter_base {

    input = [
        [ id:'test'], // meta map
        file(params.niche_test['tnfilter']['vcf'], checkIfExists: true),
        file(params.niche_test['tnfilter']['tbi'], checkIfExists: true),
        file(params.niche_test['tnfilter']['stats'], checkIfExists: true),
        [],
        [],
        []
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}

workflow test_sentieon_tnfilter_with_files {

    input = [
        [ id:'test'], // meta map
        file(params.niche_test['tnfilter']['vcf'], checkIfExists: true),
        file(params.niche_test['tnfilter']['tbi'], checkIfExists: true),
        file(params.niche_test['tnfilter']['stats'], checkIfExists: true),
        [ file(params.niche_test['tnfilter']['contamination_data'], checkIfExists: true) ],
        [ file(params.niche_test['tnfilter']['segments'], checkIfExists: true) ],
        [ file(params.niche_test['tnfilter']['orientation_data'], checkIfExists: true) ]
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}

workflow test_sentieon_tnfilter_base_stubs {

    input = [
        [ id:'test'], // meta map
        file(params.niche_test['tnfilter']['vcf'], checkIfExists: true),
        file(params.niche_test['tnfilter']['tbi'], checkIfExists: true),
        file(params.niche_test['tnfilter']['stats'], checkIfExists: true),
        [],
        [],
        []
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}
