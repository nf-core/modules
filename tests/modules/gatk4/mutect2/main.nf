#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MUTECT2 } from '../../../../modules/gatk4/mutect2/main.nf' addParams( options: [:] )

workflow test_gatk4_mutect2 {

input = [ [ id:'test', run_pon:true, run_single:false, which_norm:'1234N.recal.bam' ], // meta map
           [ file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/recalbam/1234N.recal.bam'), file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/recalbam/9876T.recal.bam') ],
           [ file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/recalbam/1234N.recal.bai'), file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/recalbam/9876T.recal.bai') ]]


    fasta = file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/reference/human_g1k_v37_decoy.small.fasta')
    fastaidx = file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/reference/human_g1k_v37_decoy.small.fasta.fai')
    dict = file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/reference/human_g1k_v37_decoy.small.fasta.dict')
    germline_resource = file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/reference/gnomAD.r2.1.1.GRCh37.small.PASS.AC.AF.only.vcf.gz')
    panel_of_normals = file('https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/reference/dbsnp_138.b37.small.vcf.gz')

    GATK4_MUTECT2 ( input )
}
