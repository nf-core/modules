#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.read_structure = "+T 12M11S+T"

include { FGBIO_FASTQTOBAM } from '../../../../modules/fgbio/fastqtobam/main.nf'

workflow test_fgbio_fastqtobam {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true), 
                file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true) ] 
              ]

    FGBIO_FASTQTOBAM ( input, "${params.read_structure}" )
}
