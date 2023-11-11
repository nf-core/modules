#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX                         } from '../../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_ALN                           } from '../../../../../modules/nf-core/bwa/aln/main.nf'
include { BWA_SAMPE                         } from '../../../../../modules/nf-core/bwa/sampe/main.nf'
include { CIRCULARMAPPER_CIRCULARGENERATOR  } from '../../../../../modules/nf-core/circularmapper/circulargenerator/main.nf'
include { CIRCULARMAPPER_REALIGNSAMFILE     } from '../../../../../modules/nf-core/circularmapper/realignsamfile/main.nf'

workflow test_circularmapper_realignsamfile {

    Channel
        .fromList(
            [
                [ id:'test', single_end:false ],
                [ [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ] ]
            ]
        ).collect()
        .set { input }
    
    fasta = [
        [ id:'test', circularextension:500, circulartarget:'MT192765.1' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    CIRCULARMAPPER_CIRCULARGENERATOR ( fasta )
    BWA_INDEX ( CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta )
    BWA_ALN ( input, BWA_INDEX.out.index )
    BWA_SAMPE ( input.join(BWA_ALN.out.sai), BWA_INDEX.out.index )
    CIRCULARMAPPER_REALIGNSAMFILE ( BWA_SAMPE.out.bam, CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta.join(CIRCULARMAPPER_CIRCULARGENERATOR.out.elongated) )
}