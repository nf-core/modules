#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_VARIANTRECALIBRATOR } from '../../../../modules/gatk4/variantrecalibrator/main.nf' addParams( options: [:] )

workflow test_gatk4_variantrecalibrator {

    input = [ [ id:'test', single_end:false ], // meta map
              file('/home/AD/gmackenz/test_storage/output/pytest_workflow_wuk87pig/gatk_tumor_only_somatic_variant_calling/output/mutect2/test.vcf', checkIfExists: true),
              file('/home/AD/gmackenz/test_storage/output/pytest_workflow_wuk87pig/gatk_tumor_only_somatic_variant_calling/output/mutect2/test.vcf.gz.tbi', checkIfExists: true)
            ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    allelespecific = false
    resource_vcfs =   [ file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
                      ]
    resource_tbis =   [ file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'], checkIfExists: true),
                        file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)
                      ]
    resource_labels = [ 'gnomAD,known=false,training=false,truth=true,prior=15.0 gnomAD.r2.1.1.vcf.gz',
                        '1000G,known=false,training=true,truth=false,prior=12.0 mills_and_1000G.indels.vcf.gz',
                        'dbsnp,known=true,training=false,truth=false,prior=10.0 dbsnp_146.hg38.vcf.gz'
                      ]
    annotation = ['DP']
    mode = 'SNP'
    create_rscript = true

    GATK4_VARIANTRECALIBRATOR ( input, fasta, fai, dict, allelespecific, resource_vcfs, resource_tbis, resource_labels, annotation, mode, create_rscript)
}
