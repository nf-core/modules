#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIMAP2_ALIGN } from '../../../../modules/minimap2/align/main.nf'

workflow test_minimap2_align_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    sam_format               = true
    preset_pacbio_reads      = false
    preset_nanopore_reads    = false
    preset_pacbio_hifi_reads = false
    preset_pacbio_overlap    = false
    preset_nanopore_overlap  = false
    preset_asm5              = false
    preset_asm10             = false
    preset_asm20             = false
    preset_nanopore_spliced  = false
    preset_pacbio_spliced    = false
    preset_short_read        = true


    MINIMAP2_ALIGN ( input, fasta, sam_format, preset_pacbio_reads, preset_nanopore_reads, preset_pacbio_hifi_reads, preset_pacbio_overlap,
                     preset_nanopore_overlap, preset_asm5, preset_asm10, preset_asm20, preset_nanopore_spliced, preset_pacbio_spliced, preset_short_read )
}

workflow test_minimap2_align_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    sam_format               = true
    preset_pacbio_reads      = false
    preset_nanopore_reads    = false
    preset_pacbio_hifi_reads = false
    preset_pacbio_overlap    = false
    preset_nanopore_overlap  = false
    preset_asm5              = false
    preset_asm10             = false
    preset_asm20             = false
    preset_nanopore_spliced  = false
    preset_pacbio_spliced    = false
    preset_short_read        = true

    MINIMAP2_ALIGN ( input, fasta, sam_format, preset_pacbio_reads, preset_nanopore_reads, preset_pacbio_hifi_reads, preset_pacbio_overlap,
                     preset_nanopore_overlap, preset_asm5, preset_asm10, preset_asm20, preset_nanopore_spliced, preset_pacbio_spliced, preset_short_read )
}
