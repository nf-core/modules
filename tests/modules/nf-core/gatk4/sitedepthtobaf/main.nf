#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SITEDEPTHTOBAF      } from '../../../../../modules/nf-core/gatk4/sitedepthtobaf/main.nf'
include { GATK4_COLLECTSVEVIDENCE   } from '../../../../../modules/nf-core/gatk4/collectsvevidence/main.nf'

workflow test_gatk4_sitedepthtobaf {

    vcf = file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)
    tbi = file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)

    input = Channel.of([
        [ id:'tumor', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        vcf,
        tbi
    ],
    [
        [ id:'normal', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        vcf,
        tbi
    ]
    )

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)


    GATK4_COLLECTSVEVIDENCE ( input, fasta, fasta_fai, dict )

    sitedepthtobaf_input = GATK4_COLLECTSVEVIDENCE.out.site_depths
                                    .combine(GATK4_COLLECTSVEVIDENCE.out.site_depths_index, by:0)
                                    .map({ meta, file, index -> [ [id:'test'], file, index ]}).groupTuple()

    GATK4_SITEDEPTHTOBAF ( sitedepthtobaf_input, [ vcf, tbi ], fasta, fasta_fai, dict )
}
