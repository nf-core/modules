#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRIDSS_GRIDSSSOMATICFILTER } from '../../../../../modules/nf-core/gridss/gridsssomaticfilter/main.nf'
include { BWA_INDEX     } from '../../../../../modules/nf-core/bwa/index/main.nf'

workflow test_gridss_gridsssomaticfilter {

    vcf = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_sv_vcf'], checkIfExists: true)
    ]
    fasta = [[id:'fasta'],file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]
    fasta_fai = [[id:'fasta_fai'],file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]

    bwa_index = BWA_INDEX(fasta).index

    GRIDSS_GRIDSSSOMATICFILTER ( vcf,[[],[]],[[],[]],[[],[]],[[],[]],[[],[]])
}
