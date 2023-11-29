#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PRINTSVEVIDENCE     } from '../../../../../modules/nf-core/gatk4/printsvevidence/main.nf'
include { GATK4_COLLECTSVEVIDENCE   } from '../../../../../modules/nf-core/gatk4/collectsvevidence/main.nf'

workflow test_gatk4_printsvevidence {

    input = Channel.of(
        [
            [ id:'normal' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ],
        [
            [ id:'tumor' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ]
    )

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)


    GATK4_COLLECTSVEVIDENCE(
        input,
        fasta,
        fasta_fai,
        dict
    )

    files = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence

    indices = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence_index

    input_svevidence = files.combine(indices, by:0).map({ meta, file, index ->
                            [ [id:'pe_files'], file, index ]
                        })
                        .groupTuple()

    GATK4_PRINTSVEVIDENCE(
        input_svevidence,
        [],
        fasta,
        fasta_fai,
        dict
    )
}

workflow test_gatk4_printsvevidence_bed_no_fasta {

    input = Channel.of(
        [
            [ id:'normal', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ],
        [
            [ id:'tumor', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
            [],
            []
        ]
    )

    fasta = []
    fasta_fai = []
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    GATK4_COLLECTSVEVIDENCE(
        input,
        fasta,
        fasta_fai,
        dict
    )

    files = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence

    indices = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence_index

    input_svevidence = files.combine(indices, by:0).map({ meta, file, index ->
                            [ [id:'pe_files'], file, index ]
                        })
                        .groupTuple()

    GATK4_PRINTSVEVIDENCE(
        input_svevidence,
        [],
        fasta,
        fasta_fai,
        dict
    )

}
