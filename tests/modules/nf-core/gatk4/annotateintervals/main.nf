#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_ANNOTATEINTERVALS } from '../../../../../modules/nf-core/gatk4/annotateintervals/main.nf'

workflow test_gatk4_annotateintervals_one_bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    dict = [[:], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_ANNOTATEINTERVALS (
        input,
        fasta,
        fasta_fai,
        dict,
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]]
    )
}

workflow test_gatk4_annotateintervals_multi_bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
        ]
    ]

    fasta = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    dict = [[:], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_ANNOTATEINTERVALS (
        input,
        fasta,
        fasta_fai,
        dict,
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]]
    )
}

workflow test_gatk4_annotateintervals_interval_list {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    fasta = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    dict = [[:], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_ANNOTATEINTERVALS (
        input,
        fasta,
        fasta_fai,
        dict,
        [[:],[]],
        [[:],[]],
        [[:],[]],
        [[:],[]]
    )
}

workflow test_gatk4_annotateintervals_mappable_regions {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    mappable_regions = [[:], file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true)]
    mappable_regions_tbi = [[:], file(params.test_data['homo_sapiens']['genome']['genome_bed_gz_tbi'], checkIfExists: true)]

    fasta = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    dict = [[:], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_ANNOTATEINTERVALS (
        input,
        fasta,
        fasta_fai,
        dict,
        mappable_regions,
        mappable_regions_tbi,
        [[:],[]],
        [[:],[]]
    )
}

workflow test_gatk4_annotateintervals_segmental_duplication_regions {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    segmental_duplication_regions = [[:], file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true)]
    segmental_duplication_regions_tbi = [[:], file(params.test_data['homo_sapiens']['genome']['genome_bed_gz_tbi'], checkIfExists: true)]

    fasta = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[:], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]
    dict = [[:], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)]

    GATK4_ANNOTATEINTERVALS (
        input,
        fasta,
        fasta_fai,
        dict,
        [[:],[]],
        [[:],[]],
        segmental_duplication_regions,
        segmental_duplication_regions_tbi
    )
}
