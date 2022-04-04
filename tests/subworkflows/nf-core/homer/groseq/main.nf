#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_GROSEQ as HOMER_GROSEQ_BAM
          HOMER_GROSEQ as HOMER_GROSEQ_BED } from '../../../../../subworkflows/nf-core/homer/groseq/main'

workflow test_homer_groseq_bam {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)]]
    def fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    HOMER_GROSEQ_BAM ( input, fasta )
}

workflow test_homer_groseq_bed {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)]]
    def fasta = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    HOMER_GROSEQ_BED ( input, fasta )
}
