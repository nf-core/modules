#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNIFFLES } from '../../../../modules/nf-core/sniffles/main.nf'

workflow test_sniffles {

    bam_tuple_ch = Channel.of([ [ id:'test' ], // meta map
                               file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true),
                               file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true),
                                ])

    fasta = Channel.of( [ [id: 'fasta'],
                         file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                        ])

    SNIFFLES ( bam_tuple_ch, fasta )
}
