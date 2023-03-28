#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_BOOT } from '../../../../../modules/nf-core/freyja/boot/main.nf'
include { FREYJA_UPDATE } from '../../../../../modules/nf-core/freyja/update/main.nf'
include { FREYJA_VARIANTS } from '../../../../../modules/nf-core/freyja/variants/main.nf'

workflow test_freyja_boot {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    FREYJA_VARIANTS (input,fasta )
    db_name= "freyja_db"
    FREYJA_UPDATE(db_name)

    variants= FREYJA_VARIANTS.out.variants
    repeats = 100
    barcodes=FREYJA_UPDATE.out.barcodes
    lineages_meta=FREYJA_UPDATE.out.lineages_meta

    FREYJA_BOOT (variants,  repeats, barcodes, lineages_meta)
}
