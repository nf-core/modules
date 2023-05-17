#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRAPHTYPER_GENOTYPE                               } from '../../../../../modules/nf-core/graphtyper/genotype/main.nf'
include { GRAPHTYPER_GENOTYPE as GRAPHTYPER_GENOTYPE_REGION } from '../../../../../modules/nf-core/graphtyper/genotype/main.nf'
include { SAMTOOLS_VIEW                                     } from '../../../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_INDEX                                    } from '../../../../../modules/nf-core/samtools/index/main.nf'

workflow test_graphtyper_genotype_single {
    
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
    ]
    reference = [ [ id: 'ref' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    ref_index = [ [ id: 'ref_index' ], file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]
    region = file(params.test_data['sarscov2']['genome']['regions_txt'], checkIfExists: true)

    GRAPHTYPER_GENOTYPE ( input, reference, ref_index, region )
}

workflow test_graphtyper_genotype_multi {
    
    input = Channel.fromList ( [
        [
            [ id: 'test1' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
        ],
        [
            [ id: 'test2' ], //meta map
            file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
        ]
    ] )
    reference = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ref_with_meta = [ [ id: 'ref' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    ref_index = [ [ id: 'ref_index' ], file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]
    region = file(params.test_data['sarscov2']['genome']['regions_txt'], checkIfExists: true)

    SAMTOOLS_VIEW ( input, reference, [])
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.cram )
    cram_grouped = SAMTOOLS_VIEW.out.cram.map{ [[id: 'group'], it[1]] }.groupTuple()
    crai_grouped = SAMTOOLS_INDEX.out.crai.map{ [[id: 'group'], it[1]] }.groupTuple()
    combined_grouped = cram_grouped.join(crai_grouped)
    GRAPHTYPER_GENOTYPE ( combined_grouped, ref_with_meta, ref_index, region )
}

workflow test_graphtyper_genotype_region {

    input = [
        [ id: 'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true) ]
    ]
    reference = [ [ id: 'ref' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)]
    ref_index = [ [ id: 'ref_index' ], file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)]
    region = []

    GRAPHTYPER_GENOTYPE_REGION ( input, reference, ref_index, region )
}

