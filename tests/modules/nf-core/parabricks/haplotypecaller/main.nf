#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_HAPLOTYPECALLER } from '../../../../../modules/nf-core/parabricks/haplotypecaller/main.nf'

workflow test_parabricks_haplotypecaller {
    input     = [
                  [ id:'input_test' ],
                  file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
                ]
    fasta     = [
                  [ id:'fasta_test' ],
                  file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
                ]

    PARABRICKS_HAPLOTYPECALLER ( input, fasta )
}

workflow test_parabricks_haplotypecaller_cram {
  input     = [
                [ id:'input_test' ],
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
              ]
  fasta     = [
                [ id:'fasta_test' ],
                file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
              ]
              
  PARABRICKS_HAPLOTYPECALLER ( input, fasta )
            
}