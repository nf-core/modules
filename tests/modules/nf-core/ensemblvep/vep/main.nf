#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_DEFAULT   } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_JSON      } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_TAB       } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF       } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_BGZIP } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_VCF_GZIP  } from '../../../../../modules/nf-core/ensemblvep/vep/main.nf'

workflow test_ensemblvep_vep_fasta_json {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_JSON ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_fasta_tab {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_TAB ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_VCF ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf_bgzip {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_VCF_BGZIP ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_fasta_vcf_gzip {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_VCF_GZIP ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    ENSEMBLVEP_VEP_DEFAULT ( input, "WBcel235", "caenorhabditis_elegans", "108", [], fasta, [] )
}

workflow test_ensemblvep_vep_no_fasta {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    ENSEMBLVEP_VEP_DEFAULT ( input, "WBcel235", "caenorhabditis_elegans", "108", [], [], [] )
}
