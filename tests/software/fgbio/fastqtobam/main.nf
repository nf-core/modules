#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.read_structure = "+T 12M11S+T"

include { FGBIO_FASTQTOBAM } from '../../../../software/fgbio/fastqtobam/main.nf' addParams( options: [args: ''] )

workflow test_fgbio_fastqtobam {

    def input = []
    input = [ [id: 'test'], //meta map
                  [ file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_1.fastq.gz'), file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/umi-dna/qiaseq/SRR7545951-small_2.fastq.gz') ] ]

    FGBIO_FASTQTOBAM ( input, "${params.read_structure}" )
}
