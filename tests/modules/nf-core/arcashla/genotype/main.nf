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
        [file("/home-link/kymmp01/hiwi/test_out/out/NA12878.mini.extracted.1.fq.gz", checkIfExists: true),
        file("/home-link/kymmp01/hiwi/test_out/out/NA12878.mini.extracted.2.fq.gz", checkIfExists: true)]
    ]

    ARCASHLA_GENOTYPE ( input )
}
