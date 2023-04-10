#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM   } from '../../../../modules/nf-core/bwa/mem/main.nf'
include { BWA_INDEX as BWA_INDEX_COV2 } from '../../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM as BWA_MEM_COV2   } from '../../../../modules/nf-core/bwa/mem/main.nf'

include { SAMTOOLS_SORT } from '../../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_COV2   } from '../../../../modules/nf-core/samtools/sort/main.nf'

include { BAMCMP } from '../../../../modules/nf-core/bamcmp/main.nf'

workflow test_bamcmp {

    input = [
        [ id:'test'], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    fasta1 = [
        [ id:'homo_sapiens_genome'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fasta2 = [
        [ id:'sarscov2_genome'], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta1 )
    BWA_MEM ( input, BWA_INDEX.out.index, false )
    SAMTOOLS_SORT (BWA_MEM.out.bam)


    BWA_INDEX_COV2 ( fasta2 )
    BWA_MEM_COV2 ( input, BWA_INDEX_COV2.out.index, false )
    SAMTOOLS_SORT_COV2 (BWA_MEM_COV2.out.bam)

    BAMCMP (SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT_COV2.out.bam, by: [0]))

}
