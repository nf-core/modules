#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_INTERSECT } from '../../../../software/bedtools/intersect/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_intersect {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
            ]
    
    output_suffix = 'bed'

    BEDTOOLS_INTERSECT ( input, output_suffix )
}

workflow test_bedtools_intersect_bam {
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    
    output_suffix = 'bam'

    BEDTOOLS_INTERSECT ( input, output_suffix )
}
