#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_VARIANTRECALIBRATOR } from '../../../../modules/gatk4/variantrecalibrator/main.nf' addParams( options: [:] )

workflow test_gatk4_variantrecalibrator {

    input = [ [ id:'test' ], // meta map
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz', checkIfExists: true),
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz.tbi', checkIfExists: true)
            ]

    fasta = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta', checkIfExists: true)
    fai = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta.fai', checkIfExists: true)
    dict = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.dict', checkIfExists: true)
    allelespecific = false
    resources = [
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz', checkIfExists: true)
                 ],
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz.tbi', checkIfExists: true)
                 ],
                 [
                 'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38_chr21.vcf.gz',
                 'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38_chr21.vcf.gz',
                 '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz',
                 'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38_chr21.vcf.gz'
                 ]
                ]
    annotation = ['QD', 'MQ', 'FS', 'SOR']
    mode = 'SNP'
    create_rscript = false

    GATK4_VARIANTRECALIBRATOR ( input, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript)
}

workflow test_gatk4_variantrecalibrator_allele_specific {

    input = [ [ id:'test' ], // meta map
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz', checkIfExists: true),
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz.tbi', checkIfExists: true)
            ]

    fasta = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta', checkIfExists: true)
    fai = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta.fai', checkIfExists: true)
    dict = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.dict', checkIfExists: true)
    allelespecific = true
    resources = [
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz', checkIfExists: true)
                 ],
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz.tbi', checkIfExists: true)
                 ],
                 [
                 'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38_chr21.vcf.gz',
                 'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38_chr21.vcf.gz',
                 '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz',
                 'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38_chr21.vcf.gz'
                 ]
                ]
    annotation = ['QD', 'MQ', 'FS']
    mode = 'SNP'
    create_rscript = false

    GATK4_VARIANTRECALIBRATOR ( input, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript)
}

workflow test_gatk4_variantrecalibrator_indel_mode {

    input = [ [ id:'test' ], // meta map
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz', checkIfExists: true),
              file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103_snpEff.ann.vcf.gz.tbi', checkIfExists: true)
            ]

    fasta = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta', checkIfExists: true)
    fai = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta.fai', checkIfExists: true)
    dict = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.dict', checkIfExists: true)
    allelespecific = false
    resources = [
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz', checkIfExists: true)
                 ],
                 [
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/hapmap_3.3.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_omni2.5.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz.tbi', checkIfExists: true),
                 file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/gatkbundle/dbsnp_138.hg38_chr21.vcf.gz.tbi', checkIfExists: true)
                 ],
                 [
                 'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38_chr21.vcf.gz',
                 'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38_chr21.vcf.gz',
                 '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38_chr21.vcf.gz',
                 'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38_chr21.vcf.gz'
                 ]
                ]
    annotation = ['QD', 'MQ', 'FS', 'SOR']
    mode = 'INDEL'
    create_rscript = false

    GATK4_VARIANTRECALIBRATOR ( input, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript)
}
