#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../../modules/nf-core/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN } from '../../../../../modules/nf-core/bowtie2/align/main.nf'
include { SEMIBIN_SINGLEEASYBIN } from '../../../../../modules/nf-core/semibin/singleeasybin/main.nf'

workflow test_semibin_singleeasybin {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    
    fasta = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
    save_unaligned = false
    sort = true

    BOWTIE2_BUILD ( fasta )
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index, save_unaligned, sort )
    
    Channel.fromPath(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
        .map { it -> [[ id:'test', single_end:true ], it] }
        .join(BOWTIE2_ALIGN.out.aligned)
        .set { input_semibin }

    SEMIBIN_SINGLEEASYBIN ( input_semibin )
   
}
