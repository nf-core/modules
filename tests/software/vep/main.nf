#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VEP } from '../../../software/vep/main.nf' addParams( vep_tag: '104.3.WBcel235', use_cache: false )

workflow test_vep {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]
            ]
    VEP ( input, "WBcel235", "caenorhabditis_elegans", "104", [] )
}
