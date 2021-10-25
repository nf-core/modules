#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZAUTILS_SEQZBINNING } from '../../../../modules/sequenzautils/seqzbinning/main.nf' addParams( options: [suffix:'.binned'] )
include { SEQUENZAUTILS_BAM2SEQZ    } from '../../../../modules/sequenzautils/bam2seqz/main.nf'    addParams( options: [:] )
include { SEQUENZAUTILS_GCWIGGLE    } from '../../../../modules/sequenzautils/gcwiggle/main.nf'    addParams( options: [:] )

workflow test_sequenzautils_seqzbinning {
    
    tumourbam = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    normalbam = file(params.test_data['homo_sapiens']['illumina']['test_paired_end_markduplicates_sorted_bam'], checkIfExists: true)
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    input_gcwiggle = [ [ id:'test' ], fasta ]
    SEQUENZAUTILS_GCWIGGLE(input_gcwiggle)

    input_bam2seqz = [ [ id:'test' ], tumourbam, normalbam ]
    SEQUENZAUTILS_BAM2SEQZ ( input_bam2seqz, fasta, SEQUENZAUTILS_GCWIGGLE.out.wig.map { it[1] } )

    SEQUENZAUTILS_SEQZBINNING ( SEQUENZAUTILS_BAM2SEQZ.out[0] )
}
