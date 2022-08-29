#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAHA_INDEX } from '../../../../modules/yaha/index/main.nf'
include { GUNZIP     } from '../../../../modules/gunzip/main.nf'
include { YAHA_ALIGN } from '../../../../modules/yaha/align/main.nf'

workflow test_yaha_align {

    input = Channel.of(
        [[ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)]
    )

    fasta = Channel.of(
        [[id: 'homo_sapien'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    )


    YAHA_INDEX ( fasta )

    GUNZIP( input )

    YAHA_ALIGN (
        GUNZIP.out.gunzip,
        YAHA_INDEX.out.nib2
            .map{ meta, nib2 -> nib2 },
        YAHA_INDEX.out.index
            .map{ meta, index -> index } )
}
