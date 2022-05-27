#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BAM_QC } from '../../../../subworkflows/local/bam_qc/main'


workflow test_bam_qc {

    bam = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'],                 checkIfExists: true)
    gff   = file(params.test_data['homo_sapiens']['genome']['genome_gff3'],                  checkIfExists: true)

    BAM_QC ( bam, fasta, gff )

}
