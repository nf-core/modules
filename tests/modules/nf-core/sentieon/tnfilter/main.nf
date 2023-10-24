#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_TNFILTER } from '../../../../../modules/nf-core/sentieon/tnfilter/main.nf'

workflow test_sentieon_tnfilter_base {

    input = [
        [ id:'test'], // meta map
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.tbi', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.stats', checkIfExists: true),
        [],
        [],
        []
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}

workflow test_sentieon_tnfilter_with_files {

    input = [
        [ id:'test'], // meta map
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.tbi', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.stats', checkIfExists: true),
        [ file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.contamination_data.tsv', checkIfExists: true) ],
        [ file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.segments', checkIfExists: true) ],
        [ file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.orientation_data.tsv', checkIfExists: true) ]
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}

workflow test_sentieon_tnfilter_base_stubs {

    input = [
        [ id:'test'], // meta map
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.tbi', checkIfExists: true),
        file('/home/ubuntu/test_data_tnhaplotyper2/sample4_vs_sample3.tnhaplotyper2.vcf.gz.stats', checkIfExists: true),
        [],
        [],
        []
    ]

    fasta = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]
    fai = [ [ id:'genome' ], // meta map
            file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    SENTIEON_TNFILTER ( input, fasta, fai )
}
