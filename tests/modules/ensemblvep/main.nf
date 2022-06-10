#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP } from '../../../modules/ensemblvep/main.nf'

workflow test_ensemblvep_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP ( input, "WBcel235", "caenorhabditis_elegans", "104", [], fasta, [] )
}

workflow test_ensemblvep_no_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    ENSEMBLVEP ( input, "WBcel235", "caenorhabditis_elegans", "104", [], [], [] )
}
