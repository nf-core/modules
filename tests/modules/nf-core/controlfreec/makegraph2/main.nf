#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONTROLFREEC_MAKEGRAPH2 } from '../../../../../modules/nf-core/controlfreec/makegraph2/main.nf'
include { CONTROLFREEC_FREEC      } from '../../../../../modules/nf-core/controlfreec/freec/main.nf'
include { UNTAR                   } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_controlfreec_makegraph2 {

    input = [
        [ id:'test', single_end:false, sex:'XX' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_mpileup'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_mpileup'], checkIfExists: true),
        [],[],[],[]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    dbsnp = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    dbsnp_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    chrfiles = [ [], file(params.test_data['homo_sapiens']['genome']['genome_21_chromosomes_dir'], checkIfExists: true) ]
    target_bed = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    UNTAR(chrfiles)
    CONTROLFREEC_FREEC (input,
                        fasta,
                        fai,
                        [],
                        dbsnp,
                        dbsnp_tbi,
                        UNTAR.out.untar.map{ it[1] },
                        [],
                        target_bed,
                        []
                        )

    makegraph_in = CONTROLFREEC_FREEC.out.ratio.join(CONTROLFREEC_FREEC.out.BAF)
    CONTROLFREEC_MAKEGRAPH2 ( makegraph_in )
}

workflow test_controlfreec_makegraph2_single {

    input = [
        [ id:'test', single_end:false, sex:'XX' ], // meta map
        [],
        file(params.test_data['homo_sapiens']['illumina']['test2_mpileup'], checkIfExists: true),
        [],[],[],[]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    dbsnp = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    dbsnp_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    chrfiles = [ [], file(params.test_data['homo_sapiens']['genome']['genome_21_chromosomes_dir'], checkIfExists: true) ]
    target_bed = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    UNTAR(chrfiles)
    CONTROLFREEC_FREEC (input,
                        fasta,
                        fai,
                        [],
                        dbsnp,
                        dbsnp_tbi,
                        UNTAR.out.untar.map{ it[1] },
                        [],
                        target_bed,
                        []
                        )

    makegraph_in = CONTROLFREEC_FREEC.out.ratio.join(CONTROLFREEC_FREEC.out.BAF)
    CONTROLFREEC_MAKEGRAPH2 ( makegraph_in )
}
