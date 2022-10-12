#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MELT } from '../../../../modules/nf-core/melt/main.nf'

workflow test_melt {
    
bam_tuple_ch = Channel.of([ [ id:'test', single_end:false ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                               []
                                ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    tranposon_file = file
    genes_file = 

    MELT ( bam_tuple_ch, fasta, fai, tranposon_file, genes_file )
}
