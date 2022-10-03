#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_STATS } from '../../../../modules/seqkit/stats/main.nf'

workflow test_seqkit_stats_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    SEQKIT_STATS ( input )
}

workflow test_seqkit_stats_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    SEQKIT_STATS ( input )
}

workflow test_seqkit_stats_nanopore {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true),
    ]

    SEQKIT_STATS ( input )
}

workflow test_seqkit_stats_genome_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ]

    SEQKIT_STATS ( input )
}

workflow test_seqkit_stats_transcriptome_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true),
    ]

    SEQKIT_STATS ( input )
}
