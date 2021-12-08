#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_joint_germline_variant_calling/main'

workflow test_gatk_joint_germline_variant_calling_skip_haplotc {
    input         = [
                  [[ id:'test' ], // meta map
                  file("https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test.g.vcf.gz",     checkIfExists: true),
                  file("https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test.g.vcf.gz.tbi", checkIfExists: true)
                  ],
                  [[ id:'test2' ], // meta map
                  file("https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test2.g.vcf.gz",     checkIfExists: true),
                  file("https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test2.g.vcf.gz.tbi", checkIfExists: true)
                  ]
                  ]
    run_haplotc  = false
    run_vqsr      = true
    fasta         = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/hg38_chr21.fasta'                         , checkIfExists: true)
    fai           = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/hg38_chr21.fasta.fai'                     , checkIfExists: true)
    dict          = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/hg38_chr21.dict'                          , checkIfExists: true)
    sites         = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/dbsnp_146.hg38_chr21.vcf.gz'     , checkIfExists: true)
    sites_tbi     = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/dbsnp_146.hg38_chr21.vcf.gz.tbi' , checkIfExists: true)
    intervals     = file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/hg38_chr21.interval_list'                                                                                                 , checkIfExists: true)
    joint_id = "test_joint"
    allelespecific = false
    resources = [
                 [
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/hapmap_3.3.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/1000G_omni2.5.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/dbsnp_138.hg38_chr21.vcf.gz', checkIfExists: true)
                 ],
                 [
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/hapmap_3.3.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/1000G_omni2.5.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/dbsnp_138.hg38_chr21.vcf.gz.tbi', checkIfExists: true)
                 ],
                 [
                 'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38_chr21.vcf.gz',
                 'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38_chr21.vcf.gz',
                 '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz',
                 'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38_chr21.vcf.gz'
                 ]
                ]
    annotation = ['QD', 'FS', 'SOR']
    mode = 'SNP'
    create_rscript = false
    truthsensitivity = '99.0'
    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, run_haplotc, run_vqsr, fasta, fai, dict, sites, sites_tbi, intervals, joint_id, allelespecific, resources, annotation, mode, create_rscript, truthsensitivity )
}

workflow test_gatk_joint_germline_variant_calling_skip_vqsr {
    input         = [
                  [[ id:'test' ], // meta map
                  file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam',     checkIfExists: true),
                  file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam.bai', checkIfExists: true)
                  ],
                  [[ id:'test2' ], // meta map
                  file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam',     checkIfExists: true),
                  file('https://raw.githubusercontent.com/GCJMackenzie/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam.bai', checkIfExists: true)
                  ]
                  ]
    run_haplotc  = true
    run_vqsr      = false
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    sites = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    sites_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)
    intervals = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    joint_id = "test_joint"
    allelespecific = []
    resources = []
    annotation = []
    mode = []
    create_rscript = []
    truthsensitivity = []
    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, run_haplotc, run_vqsr, fasta, fai, dict, sites, sites_tbi, intervals, joint_id, allelespecific, resources, annotation, mode, create_rscript, truthsensitivity )
}
