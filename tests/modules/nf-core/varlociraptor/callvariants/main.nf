#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL } from '../../../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main.nf'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES as VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR } from '../../../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main.nf'
include { VARLOCIRAPTOR_PREPROCESS  as      VARLOCIRAPTOR_PREPROCESS_NORMAL          } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'
include { VARLOCIRAPTOR_PREPROCESS  as      VARLOCIRAPTOR_PREPROCESS_TUMOR         } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'

include { VARLOCIRAPTOR_CALLVARIANTS                } from '../../../../../modules/nf-core/varlociraptor/callvariants/main.nf'

workflow test_varlociraptor_callvariants {

    bam_normal = [
        [ id:'test_normal', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
    ]

    bam_tumor = [
        [ id:'test_tumor', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai= [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL( bam_normal, fasta, fai)
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR( bam_tumor, fasta, fai)

    input_normal = Channel.of([
        [ id:'test_normal', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
    ]).collect()

    input_tumor = Channel.of([
        [ id:'test_tumor', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_genome_vcf'], checkIfExists: true),
    ]).collect()

    input_normal.view()
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json.view()

    VARLOCIRAPTOR_PREPROCESS_NORMAL(input_normal.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_NORMAL.out.alignment_properties_json), fasta, fai)
    VARLOCIRAPTOR_PREPROCESS_TUMOR(input_tumor.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES_TUMOR.out.alignment_properties_json), fasta, fai)

    //scenario = file("./test.yml", checkIfExists: true)
    VARLOCIRAPTOR_CALLVARIANTS ( VARLOCIRAPTOR_PREPROCESS_NORMAL.out.vcf_gz, VARLOCIRAPTOR_PREPROCESS_TUMOR.out.vcf_gz, [] )
}
