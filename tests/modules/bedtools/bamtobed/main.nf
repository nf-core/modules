#!/usr/bin/env nextflow



include { BEDTOOLS_BAMTOBED } from '../../../../modules/bedtools/bamtobed/main.nf'

workflow test_bedtools_bamtobed {
    input = [ [ id:'test'], //meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
            ]

    BEDTOOLS_BAMTOBED ( input )
}
