#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEKA               } from '../../../../modules/nf-core/peka/main.nf'

workflow test_peka {

    bed_crosslinks = [ [ id:'test' ], file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/peka/chr21_HepG2-PCBP1-merged.xl.bed", checkIfExists: true) ]
    bed_peaks      = [ [ id:'test' ], file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/peka/chr21_HepG2-PCBP1-merged.xl10_200_density2_peaks.bed", checkIfExists: true) ]
    regions        = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/peka/chr21_gencode_regions.gtf", checkIfExists: true)
    fasta          = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/peka/chr21.GRCh38.p12.genome.masked.fa", checkIfExists: true)
    fai            = file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/peka/chr21.GRCh38.p12.genome.masked.fa.fai", checkIfExists: true)

    PEKA(
        bed_peaks,
        bed_crosslinks,
        fasta,
        fai,
        regions
    )
}
