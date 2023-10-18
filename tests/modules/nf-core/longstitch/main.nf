#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGSTITCH } from '../../../../modules/nf-core/longstitch/main.nf'

workflow test_longstitch {
    
    input_assembly = [
        [ id:'test_assembly', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    input_reads = [
        [ id:'test_reads', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    mode = "ntLink-arks"
    genome_size = "3e9"

    LONGSTITCH ( input_assembly, input_reads, mode, genome_size )
}
