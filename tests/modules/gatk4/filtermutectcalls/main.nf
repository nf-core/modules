#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERMUTECTCALLS } from '../../../../modules/gatk4/filtermutectcalls/main.nf' addParams( options: [:] )

workflow test_gatk4_filtermutectcalls_base {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.tbi", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.stats", checkIfExists: true),
              [],
              [],
              [],
              []]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input , fasta , fastaidx , dict )
}

workflow test_gatk4_filtermutectcalls_with_files {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.tbi", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.stats", checkIfExists: true),
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_test2_paired_mutect2_calls.artifact-prior.tar.gz", checkIfExists: true)],
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_test2_paired.segmentation.table", checkIfExists: true)],
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_test2_paired.contamination.table", checkIfExists: true)],
              []]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input , fasta , fastaidx , dict )
}

workflow test_gatk4_filtermutectcalls_use_val {

    input = [ [ id:'test'], // meta map
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.tbi", checkIfExists: true),
              file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz.stats", checkIfExists: true),
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_test2_paired_mutect2_calls.artifact-prior.tar.gz", checkIfExists: true)],
              [file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/test_test2_paired.segmentation.table", checkIfExists: true)],
              [],
              '20.0']

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fastaidx = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_FILTERMUTECTCALLS ( input , fasta , fastaidx , dict )
}
