#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE }       from '../../../../software/rsem/preparereference/main.nf'   addParams(options: [:])
include { RSEM_CALCULATEEXPRESSION }    from '../../../../software/rsem/calculateexpression/main.nf'   addParams(options: [args: "--star"])

workflow test_rsem_calculateexpression {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    input = [ [ id:'test', single_end:true ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]

    RSEM_PREPAREREFERENCE ( fasta, gtf )

    RSEM_CALCULATEEXPRESSION( input, RSEM_PREPAREREFERENCE.out.index )
}
