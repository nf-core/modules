#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CANU } from '../../../../modules/nf-core/canu/main.nf'

workflow test_hicanu {
    
    input = [
        [ id:'test', single_end:true ], // meta map
	file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    mode = "-pacbio-hifi"
    genome_size = '1k' 

    CANU ( input, mode, genome_size )
}
