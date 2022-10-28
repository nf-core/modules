#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VGAN_HAPLOCART } from '../../../../../modules/nf-core/vgan/haplocart/main.nf' addParams( options: [args: ''] )

workflow test_vgan_haplocart_paired_end_interleaved {

    input = [
        [ id:'test', single_end:false, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        []
        ]

    VGAN_HAPLOCART(input)
}

workflow test_vgan_haplocart_paired_end_separate {

    rCRS_reads = file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true)
    rCRS_reads.copyTo(workDir + 'reads2.fq.gz')
    input = [
        [ id:'test', single_end:false, format:'fastq'], // meta map
        rCRS_reads,
        file(workDir + "reads2.fq.gz", checkIfExists: true)
            ]

    VGAN_HAPLOCART(input)
}

workflow test_vgan_haplocart_single_end {
    
    input = [
        [ id:'test', single_end:true, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        []
        ]

    VGAN_HAPLOCART(input)
}


