#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_REF } from '../../../../modules/kallistobustools/ref/main.nf'

workflow test_kallistobustools_ref_standard {

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    workflow    = "standard"

    KALLISTOBUSTOOLS_REF(fasta, gtf, workflow)
}

workflow test_kallistobustools_ref_lamanno {

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    workflow    = "lamanno"

    KALLISTOBUSTOOLS_REF( fasta, gtf, workflow)
}

workflow test_kallistobustools_ref_nucleus {

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf         = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    workflow    = "nucleus"

    KALLISTOBUSTOOLS_REF( fasta, gtf, workflow)
}

