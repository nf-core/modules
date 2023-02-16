#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_FQ2BAM } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { PARABRICKS_FQ2BAM as PARABRICKS_FQ2BAM_MKDUP_NOQC } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { PARABRICKS_FQ2BAM as PARABRICKS_FQ2BAM_NOMKDUP_NOQC } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { PARABRICKS_FQ2BAM as PARABRICKS_FQ2BAM_MKDUP_QC } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { PARABRICKS_FQ2BAM as PARABRICKS_FQ2BAM_NOMKDUP_QC } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { BWA_INDEX } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_parabricks_fq2bam_pe_default {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[] )

}

workflow test_parabricks_fq2bam_se_default {
    
    input = [
        [ id:'test', single_end:true],
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)],
        []
    ]

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=[] )
}

workflow test_parabricks_fq2bam_pe_mkdup_noqc {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM_MKDUP_NOQC ( input, fasta, BWA_INDEX.out.index, known_sites=[] )
}

workflow test_parabricks_fq2bam_pe_nomkdup_noqc {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM_NOMKDUP_NOQC ( input, fasta, BWA_INDEX.out.index, known_sites=[] )
}

workflow test_parabricks_fq2bam_pe_mkdup_qc {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM_MKDUP_QC ( input, fasta, BWA_INDEX.out.index, known_sites=[] )
}

workflow test_parabricks_fq2bam_pe_nomkdup_qc {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM_NOMKDUP_QC ( input, fasta, BWA_INDEX.out.index, known_sites=[] )
}

workflow test_parabricks_fq2bam_pe_bqsr {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
    ]
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    known_sites = [
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites )
}

workflow test_parabricks_fq2bam_pe_bqsr2 {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        []
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
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites )
}


workflow test_parabricks_fq2bam_pe_bqsr_intervals {
    
    input = [
        [ id:'test', single_end:false],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ],
        [
        file(params.test_data['sarscov2']['genome']['baits_interval_list'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)
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
    PARABRICKS_FQ2BAM ( input, fasta, BWA_INDEX.out.index, known_sites=known_sites )
}