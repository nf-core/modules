#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_FQ2BAM } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { BWA_INDEX } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_parabricks_fq2bam_pe_mkdup {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[], interval_file=[], markdups=true, qc_metrics=false )

}

workflow test_parabricks_fq2bam_pe_nomkdup {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[], interval_file=[], markdups=false, qc_metrics=false )
}

workflow test_parabricks_fq2bam_se_mkdup {
    
    input = [
        [ id:'test', single_end:true],
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[], interval_file=[], markdups=true, qc_metrics=false )
}

workflow test_parabricks_fq2bam_se_mkdup_qc {
    
    input = [
        [ id:'test', single_end:true],
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[], interval_file=[], markdups=true, qc_metrics=true )
}


workflow test_parabricks_fq2bam_pe_mkdup_bqsr {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    known_sites = [
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites, interval_file=[], markdups=true, qc_metrics=false )
}

workflow test_parabricks_fq2bam_pe_mkdup_bqsr2 {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    known_sites = [
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites, interval_file=[], markdups=true, qc_metrics=false )
}


workflow test_parabricks_fq2bam_pe_mkdup_bqsr_intervals {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    known_sites = [
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ]

    interval_file = [
        file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites, interval_file=interval_file, markdups=true, qc_metrics=false )
}