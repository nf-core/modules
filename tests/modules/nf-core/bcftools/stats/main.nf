#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../../../../../modules/nf-core/bcftools/stats/main.nf'

workflow test_bcftools_stats {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions   = [ [], [] ]
    targets   = [ [], [] ]
    samples   = [ [], [] ]
    exons     = [ [], [] ]
    reference = [ [], [] ]

    BCFTOOLS_STATS ( input, regions, targets, samples, exons, reference)
}

workflow test_bcftools_stats_regions {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]
    regions   = [ [id:'test'], file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true) ]
    targets   = [ [], [] ]
    samples   = [ [], [] ]
    exons     = [ [], [] ]
    reference = [ [], [] ]

    BCFTOOLS_STATS ( input, regions, targets, samples, exons, reference)
}

workflow test_bcftools_stats_targets {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions   = [ [], [] ]
    targets   = [ [id:'test'], file(params.test_data['sarscov2']['illumina']['test2_vcf_targets_tsv_gz'], checkIfExists: true) ]
    samples   = [ [], [] ]
    exons     = [ [], [] ]
    reference = [ [], [] ]

    BCFTOOLS_STATS ( input, regions, targets, samples, exons, reference)
}

workflow test_bcftools_stats_exons {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions   = [ [], [] ]
    targets   = [ [], [] ]
    samples   = [ [], [] ]
    exons     = [ [id: "exon_test"], file("/Users/lamnidis/Software/github/TCLamnidis/test-datasets/data/delete_me/bcftools/stats/exons.tsv.gz") ]
    reference = [ [], [] ]

    BCFTOOLS_STATS ( input, regions, targets, samples, exons, reference)
}

workflow test_bcftools_stats_reference {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions   = [ [], [] ]
    targets   = [ [], [] ]
    samples   = [ [], [] ]
    exons     = [ [], [] ]
    reference = [ [id: 'ref_test'], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BCFTOOLS_STATS ( input, regions, targets, samples, exons, reference)
}
