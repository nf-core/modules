#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_VARIANT_CALLING_FREEBAYES } from '../../../../subworkflows/nf-core/bam_variant_calling_freebayes/main.nf'

workflow test_bam_variant_calling_freebayes {
    
    input = [ [ id: "test" ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
              [],
              [],
              []
            ]

    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)


    BAM_VARIANT_CALLING_FREEBAYES ( input, genome, [], [], [] )
}
