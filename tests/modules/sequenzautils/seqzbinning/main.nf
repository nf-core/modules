#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZAUTILS_SEQZBINNING } from '../../../../modules/sequenzautils/seqzbinning/main.nf' addParams( options: [:] )

workflow test_sequenzautils_seqzbinning {
    
    tumourbam = file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    normalbam = file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)

    input = [ [ id:'test' ], // meta map
              tumourbam,
              normalbam
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    wig   = file(params.test_data['sarscov2']['illumina']['test_wig_gz'], checkIfExists: true)

    SEQUENZAUTILS_BAM2SEQZ ( input, fasta, wig )

    SEQUENZAUTILS_SEQZBINNING ( SEQUENZAUTILS_BAM2SEQZ.out[0] )
}
