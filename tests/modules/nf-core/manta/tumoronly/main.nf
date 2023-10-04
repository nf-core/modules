#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MANTA_TUMORONLY } from '../../../../../modules/nf-core/manta/tumoronly/main.nf'

workflow test_manta_tumoronly {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        [], []
    ]

    fasta = [ [id:'genome'],
              file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
            ]

    fai   = [ [id:'genome'],
              file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
            ]

    config = Channel.of("[manta]", "enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0")
        .collectFile(name:"manta_options.ini", newLine:true)

    MANTA_TUMORONLY ( input, fasta, fai, config )
}

workflow test_manta_tumoronly_target_bed {
    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed_gz_tbi'], checkIfExists: true)
    ]

    fasta = [ [id:'genome'],
              file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
            ]

    fai   = [ [id:'genome'],
              file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
            ]

    MANTA_TUMORONLY ( input, fasta, fai, [] )
}
