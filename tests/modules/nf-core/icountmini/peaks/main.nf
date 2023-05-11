#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ICOUNTMINI_SEGMENT } from '../../../../../modules/nf-core/icountmini/segment/main.nf'
include { ICOUNTMINI_SIGXLS } from '../../../../../modules/nf-core/icountmini/sigxls/main.nf'
include { ICOUNTMINI_PEAKS  } from '../../../../../modules/nf-core/icountmini/peaks/main.nf'

workflow test_icountmini_peaks {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_gencode_gtf'], checkIfExists: true)
    ]

    ICOUNTMINI_SEGMENT ( 
        input,
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    )

    bed = [
        [  id:'test' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/crosslinks/clippy.bed", checkIfExists: true)
    ]
    gtf = ICOUNTMINI_SEGMENT.out.gtf.flatten().last()

    ICOUNTMINI_SIGXLS (
        bed,
        gtf
    )

    peaks_input = bed
        .map{ [ it[0].id, it[0], it[1] ] }
        .join( ICOUNTMINI_SIGXLS.out.sigxls.map{ [ it[0].id, it[0], it[1] ] } )
        .map { [ it[1], it[2], it[4] ] }

    ICOUNTMINI_PEAKS (
        peaks_input
    )
}
