#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF   } from '../../../../../modules/nf-core/kallistobustools/ref/main.nf'
include { KALLISTOBUSTOOLS_COUNT } from '../../../../../modules/nf-core/kallistobustools/count/main.nf'

workflow test_kallistobustools_count {

    input   = [
        [id:'test'], // meta map
        [ 
          file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_gex1_fastq_1_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_10k_pbmc_cmo_gex1_fastq_2_gz'], checkIfExists: true)
        ]
    ]

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    sc_workflow = "standard"
    technology  = "10XV3"

    KALLISTOBUSTOOLS_REF(fasta, gtf, sc_workflow)
    KALLISTOBUSTOOLS_COUNT ( 
      input, 
      KALLISTOBUSTOOLS_REF.out.index, 
      KALLISTOBUSTOOLS_REF.out.t2g, 
      KALLISTOBUSTOOLS_REF.out.cdna_t2c.ifEmpty{ [] },   // when empty the module doesn't run unless something is passed. 
      KALLISTOBUSTOOLS_REF.out.intron_t2c.ifEmpty{ [] }, // when empty the module doesn't run unless something is passed.
      technology 
    )
}
