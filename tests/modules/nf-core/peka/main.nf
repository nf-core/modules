#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEKA               } from '../../../../modules/nf-core/peka/main.nf'
include { CLIPPY             } from '../../../../modules/nf-core/clippy/main.nf'
include { ICOUNTMINI_SEGMENT } from '../../../../modules/nf-core/icountmini/segment/main.nf'
include { GUNZIP             } from '../../../../modules/nf-core/gunzip/main.nf'

workflow test_peka {
    input = [
        [  id:'test' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/crosslinks/clippy.bed", checkIfExists: true)
    ]

    CLIPPY ( 
        input,
        file(params.test_data['homo_sapiens']['genome']['genome_21_gencode_gtf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
        )

    seg_input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_gencode_gtf'], checkIfExists: true)
    ]

    ICOUNTMINI_SEGMENT ( 
        seg_input,
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
        )

    bed_crosslinks = [ [  id:'test' ], file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/crosslinks/clippy.bed", checkIfExists: true) ]
    bed_peaks      = CLIPPY.out.peaks
    fasta          = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai            = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    GUNZIP ( ICOUNTMINI_SEGMENT.out.regions )

    PEKA(
        bed_peaks,
        bed_crosslinks,
        GUNZIP.out.gunzip.map{ it[1] },
        fasta,
        fai
    )
}
