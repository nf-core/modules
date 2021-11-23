#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYVQSR } from '../../../../modules/gatk4/applyvqsr/main.nf' addParams( options: [:] )

workflow test_gatk4_applyvqsr {
    input = [ [ id:'test'], // meta map
              [ file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103.vcf.gz', checkIfExists: true)],
              [ file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/germline_calling/variants/HaplotypeCaller_disease_103.vcf.gz.tbi', checkIfExists: true)],
              file('/home/AD/gmackenz/test_data/francesco/variantrecals/base/test.recal', checkIfExists: true),
              file('/home/AD/gmackenz/test_data/francesco/variantrecals/base/test.recal.idx', checkIfExists: true),
              file('/home/AD/gmackenz/test_data/francesco/variantrecals/base/test.tranches', checkIfExists: true)
            ]
    fasta = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta', checkIfExists: true)
    fai = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.fasta.fai', checkIfExists: true)
    dict = file('https://raw.githubusercontent.com/lescai-teaching/datasets_class/master/reference/sequence/Homo_sapiens_assembly38_chr21.dict', checkIfExists: true)
    allelespecific = false
    truthsensitivity = '99'
    mode = 'SNP'

    GATK4_APPLYVQSR ( input, fasta, fai, dict, allelespecific, truthsensitivity, mode )
}
