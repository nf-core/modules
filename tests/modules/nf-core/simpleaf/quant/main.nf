#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SIMPLEAF_INDEX } from '../../../../../modules/nf-core/simpleaf/index/main.nf'
include { SIMPLEAF_QUANT } from '../../../../../modules/nf-core/simpleaf/quant/main.nf'

workflow test_simpleaf_quant {
    
    genome_fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    meta = [ 'id': 'human'] 

    SIMPLEAF_INDEX ( 
        [ meta, genome_fasta ],
        [ meta, gtf ],
        [[],[]]
    )

    input = [ 
        [ id:'test_10x', single_end:false, strandedness:'auto'  ], // meta map
        '10xv3', // Chemistry
        [ 
            file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_1_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_5k_cmvpos_tcells_gex1_fastq_2_gz'], checkIfExists: true)
        ]
    ]

    SIMPLEAF_QUANT ( 
        input,
        SIMPLEAF_INDEX.out.index,
        SIMPLEAF_INDEX.out.transcript_tsv,
        'cr-like',
        [[],[]]
    )
}

