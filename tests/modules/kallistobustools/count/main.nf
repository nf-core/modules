#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_COUNT } from '../../../../modules/kallistobustools/count/main.nf'

workflow test_kallistobustools_count {

    input   = [
        [id:'test'], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    index      = file("https://github.com/FloWuenne/test-datasets/blob/scrnaseq/reference/kallistobustools/kb_ref.idx?raw=true", checkIfExists: true)
    t2g        = file("https://raw.githubusercontent.com/FloWuenne/test-datasets/scrnaseq/reference/kallistobustools/t2g.txt", checkIfExists: true)
    t1c        = []
    t2c        = []
    workflow   = "standard"
    technology = "10XV3"

    KALLISTOBUSTOOLS_COUNT ( input, index, t2g, t1c, t2c, workflow, technology )
}
