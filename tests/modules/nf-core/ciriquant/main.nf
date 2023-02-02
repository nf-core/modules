#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/nf-core/bwa/index/main.nf'
include { CIRIQUANT } from '../../../../modules/nf-core/ciriquant/main.nf'
include { HISAT2_BUILD } from '../../../modules/nf-core/hisat2/build/main.nf'
include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/nf-core/hisat2/extractsplicesites/main.nf'

workflow test_ciriquant {

    reads = [
        [ id:'test', single_end:false ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_rnaseq_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    BWA_INDEX( [[id:'test'], fasta] )
    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    CIRIQUANT( reads, gtf, fasta, BWA_INDEX.out.index.map{ meta, index -> [ index ] }, HISAT2_BUILD.out.index )

}
