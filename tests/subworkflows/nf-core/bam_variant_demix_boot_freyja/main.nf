#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_VARIANT_DEMIX_BOOT_FREYJA } from '../../../../subworkflows/nf-core/bam_variant_demix_boot_freyja/main.nf'
include { FREYJA_UPDATE                 } from '../../../../modules/nf-core/freyja/update/main.nf'

workflow test_bam_variant_demix_boot_freyja_nodb {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    repeats= 100
    barcodes= []
    lineage_meta= []
    db_name="freyja_db"
    BAM_VARIANT_DEMIX_BOOT_FREYJA ( input, fasta,repeats, db_name, barcodes, lineage_meta )
}

workflow test_bam_variant_demix_boot_freyja_withdb {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    db_name="freyja_db"
    FREYJA_UPDATE(db_name)

    repeats=100
    barcodes= FREYJA_UPDATE.out.barcodes
    lineage_meta= FREYJA_UPDATE.out.lineages_meta

    BAM_VARIANT_DEMIX_BOOT_FREYJA ( input, fasta,repeats,db_name, barcodes, lineage_meta )
}
