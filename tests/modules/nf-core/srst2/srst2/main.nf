#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRST2_SRST2 } from '../../../../../modules/nf-core/srst2/srst2/main.nf'

workflow test_srst2_srst2_exit {

    input = [
        [ id:'test', single_end:false, db:"test"], // meta map
        [ file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true) ],
//        [("")]
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/resFinder_20180221_srst2.fasta')
    ]

    SRST2_SRST2(input)
}

workflow test_srst2_srst2_mlst {

    input = [
        [ id:'test', single_end:false, db:"mlst"], // meta map
        [ file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/SRR9067271_1.fastq.gz", checkIfExists: true),
          file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/SRR9067271_2.fastq.gz", checkIfExists: true) ],
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/MLST_DB.fas')
    ]

    SRST2_SRST2(input)
}

workflow test_srst2_srst2_paired_end {

    input = [
        [ id:'test', single_end:false, db:"gene"], // meta map
        [ file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true) ],
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/resFinder_20180221_srst2.fasta')  // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    SRST2_SRST2(input)
}

workflow test_srst2_srst2_single_end {

    input = [
        [ id:'test', single_end:true, db:"gene" ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/srst2/resFinder_20180221_srst2.fasta')  // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    SRST2_SRST2(input)
}
