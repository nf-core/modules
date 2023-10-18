#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARSOLO                          } from '../../../../../modules/nf-core/star/starsolo/main.nf'
include { STAR_GENOMEGENERATE               } from '../../../../../modules/nf-core/star/genomegenerate/main.nf'

workflow test_starsolo {
    
    input = [ [ id:'test_starsolo', umi_len:'12'  ],
            'CB_UMI_Simple',
            [ file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_1_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_2_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'test_fasta' ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ]
    gtf = [
        [ id:'test_gtf' ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true) ]
    ]

    STAR_GENOMEGENERATE ( fasta, gtf )
    STARSOLO ( input, STAR_GENOMEGENERATE.out.index )
}
