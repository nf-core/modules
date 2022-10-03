#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTK_STANDARDIZE } from '../../../../../modules/nf-core/svtk/standardize/main.nf'
include { MANTA_GERMLINE   } from '../../../../../modules/nf-core/manta/germline/main.nf'

workflow test_svtk_standardize {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed_gz_tbi'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    SVTK_STANDARDIZE ( 
        MANTA_GERMLINE.out.diploid_sv_vcf,
        fasta_fai
    )
}

workflow test_svtk_standardize_no_contigs {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed_gz_tbi'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    SVTK_STANDARDIZE ( 
        MANTA_GERMLINE.out.diploid_sv_vcf,
        []
    )
}
