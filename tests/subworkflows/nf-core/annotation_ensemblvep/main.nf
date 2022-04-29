#!/usr/bin/env nextflow



include { ANNOTATION_ENSEMBLVEP } from '../../../../subworkflows/nf-core/annotation_ensemblvep/main'

workflow annotation_ensemblvep {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    ANNOTATION_ENSEMBLVEP ( input, "WBcel235", "caenorhabditis_elegans", "104", [] )
}
