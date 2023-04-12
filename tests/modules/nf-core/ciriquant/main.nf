#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/nf-core/bwa/index/main.nf'
include { CIRIQUANT } from '../../../../modules/nf-core/ciriquant/main.nf'
include { HISAT2_BUILD } from '../../../modules/nf-core/hisat2/build/main.nf'
include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/nf-core/hisat2/extractsplicesites/main.nf'

workflow test_ciriquant {

    reads = [
        [ id:'fust1', single_end:false ],
        [
            file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/fastq/fust1_rep1_1.fastq.gz"),
            file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/fastq/fust1_rep1_2.fastq.gz")
        ]
    ]

    fasta = file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.fa")
    gtf   = file("https://raw.githubusercontent.com/nf-core/test-datasets/circrna/reference/chrI.gtf")

    BWA_INDEX( [[id:'test'], fasta] )
    HISAT2_EXTRACTSPLICESITES ( gtf )
    HISAT2_BUILD( fasta, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    CIRIQUANT( reads, gtf, fasta, BWA_INDEX.out.index.map{ meta, index -> [ index ] }, HISAT2_BUILD.out.index )

}
