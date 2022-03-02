#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_INDEXFEATUREFILE } from '../../../../modules/gatk4/indexfeaturefile/main.nf'

workflow test_gatk4_indexfeaturefile_bed {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true) 
    ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_bed_gz {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed_gz'], checkIfExists: true)
    ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_vcf {
    
    input = [ 
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true) 
    ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_vcf_gz {
    
    input = [ 
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true) 
    ]

    GATK4_INDEXFEATUREFILE ( input )
}
