#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MPILEUP } from '../../../../modules/samtools/mpileup/main.nf'
include { SEQUENCETOOLS_PILEUPCALLER } from '../../../../modules/sequencetools/pileupcaller/main.nf'

workflow test_sequencetools_pileupcaller {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                []
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    snpfile = file('/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/chr_21.snp')

    SAMTOOLS_MPILEUP ( input, fasta )

    input2 = SAMTOOLS_MPILEUP.out.mpileup.combine(snpfile)
    SEQUENCETOOLS_PILEUPCALLER ( input2, 'randomHaploid', 'EIGENSTRAT' )
}
