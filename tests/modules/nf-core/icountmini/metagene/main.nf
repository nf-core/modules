#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ICOUNTMINI_SEGMENT } from '../../../../../modules/nf-core/icountmini/segment/main.nf'
include { ICOUNTMINI_METAGENE} from '../../../../../modules/nf-core/icountmini/metagene/main.nf'

workflow test_icountmini_metagene {
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

    ICOUNTMINI_METAGENE (
        bed,
        ICOUNTMINI_SEGMENT.out.gtf.map{ it[1] }
    )
}

