#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_FQ2BAM } from '../../../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { BWA_INDEX } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_parabricks_fq2bam_pe_mkdup {
    
    input = [
        [ id:'test', single_end:false, markdups:true, qc_metrics:false],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_2.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]

    PARABRICKS_FQ2BAM ( input )
}

workflow test_parabricks_fq2bam_pe_nomkdup {
    
    input = [
        [ id:'test', single_end:false, markdups:false, qc_metrics:false],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_2.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]

    PARABRICKS_FQ2BAM ( input, known_sites=[], interval_file=[], markdups=false, qc_metrics=false )
}

workflow test_parabricks_fq2bam_se {
    
    input = [
        [ id:'test', single_end:true, markdups:true, qc_metrics:false],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]

    PARABRICKS_FQ2BAM ( input, known_sites=[], interval_file=[], markdups=true, qc_metrics=false )
}


workflow test_parabricks_fq2bam_pe_mkdup_bqsr {
    
    input = [
        [ id:'test', single_end:false],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_2.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]
    known_sites = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.known_indels.vcf.gz", checkIfExists: true)
    ]

    PARABRICKS_FQ2BAM ( input, known_sites=known_sites, interval_file=[], markdups=true, qc_metrics=false )
}

workflow test_parabricks_fq2bam_pe_mkdup_bqsr2 {
    
    input = [
        [ id:'test', single_end:false],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_2.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]
    known_sites = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.known_indels.vcf.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.known_indels_2.vcf.gz", checkIfExists: true)
    ]

    PARABRICKS_FQ2BAM ( input, known_sites=known_sites, interval_file=[], markdups=true, qc_metrics=false )
}


workflow test_parabricks_fq2bam_pe_mkdup_bqsr_intervals {
    
    input = [
        [ id:'test', single_end:false ],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_1.fq.gz", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Data/sample_2M_2.fq.gz", checkIfExists: true)],
        [file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.amb", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.ann", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.bwt", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.pac", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.fasta.sa", checkIfExists: true)]
    ]
    known_sites = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/Homo_sapiens_assembly38.known_indels.vcf.gz", checkIfExists: true)
    ]
    interval_file = [
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_1.bed", checkIfExists: true),
        file("/home/bsiranos/parabricks_demo/parabricks_sample/Ref/chr1/intervals_2.bed", checkIfExists: true)
    ]
    markdups = true
    qc_metrics = true 
    PARABRICKS_FQ2BAM ( input, known_sites=known_sites, interval_file=interval_file, markdups=true, qc_metrics=true )
}