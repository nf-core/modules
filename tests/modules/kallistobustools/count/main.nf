#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF   } from '../../../../modules/kallistobustools/ref/main.nf'
include { KALLISTOBUSTOOLS_COUNT } from '../../../../modules/kallistobustools/count/main.nf'

workflow test_kallistobustools_count {

    input   = [
        [id:'test'], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_10x_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_10x_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    sc_workflow    = "standard"
    technology = "10XV3"
    use_t1c    = true
    use_t2c    = true

    KALLISTOBUSTOOLS_REF(fasta, gtf, sc_workflow)
    KALLISTOBUSTOOLS_COUNT ( input, KALLISTOBUSTOOLS_REF.out.index, KALLISTOBUSTOOLS_REF.out.t2g, KALLISTOBUSTOOLS_REF.out.cdna_t2c, KALLISTOBUSTOOLS_REF.out.intron_t2c, use_t1c, use_t2c, sc_workflow, technology )
}
