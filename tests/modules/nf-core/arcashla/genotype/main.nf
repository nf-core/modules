#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARCASHLA_GENOTYPE } from '../../../../../modules/nf-core/arcashla/genotype/main.nf'

workflow test_arcashla_genotype {
    
    // input = [
    //     [ id:'test', single_end:false ], // meta map
    //     file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    // ]
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [file("/Users/students/Documents/hiwi/hlatyping/dev_tests/test_dir/test_out/arcashla/NA12878.mini.extracted.1.fq.gz", checkIfExists: true),
        file("/Users/students/Documents/hiwi/hlatyping/dev_tests/test_dir/test_out/arcashla/NA12878.mini.extracted.2.fq.gz", checkIfExists: true)]
    ]

    ARCASHLA_GENOTYPE ( input )
}
