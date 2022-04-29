#!/usr/bin/env nextflow



include { ENSEMBLVEP } from '../../../modules/ensemblvep/main.nf'

workflow test_ensemblvep {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    ENSEMBLVEP ( input, "WBcel235", "caenorhabditis_elegans", "104", [] )
}
