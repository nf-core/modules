#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF2MAF } from '../../../../modules/nf-core/vcf2maf/main.nf'
include { UNTAR   } from '../../../../modules/nf-core/untar/main.nf'

workflow test_vcf2maf_no_vep {

    input_vcf = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]
    fasta = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]

    VCF2MAF ( input_vcf, fasta, [] )
}

workflow test_vcf2maf_vep {

    input_vcf = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)
    ]
    fasta     = [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    vep_cache = [ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['vep_cache'], checkIfExists: true) ]

    vep_cache_unpacked = UNTAR(vep_cache).untar.map { it[1] }
    VCF2MAF ( input_vcf, fasta, vep_cache_unpacked)
}
