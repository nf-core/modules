#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VGAN_HAPLOCART } from '../../../../../modules/nf-core/vgan/haplocart/main.nf' addParams( options: [args: ''] )

workflow test_vgan_haplocart_single_end {
    
    input = [
        [ id:'test', single_end:true, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        ]

    VGAN_HAPLOCART(input)
}

workflow test_vgan_haplocart_paired_end_interleaved {

    input = [
        [ id:'test', single_end:false, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        ]

    VGAN_HAPLOCART(input)
}
