#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_INTERVALLISTTOOLS AS GATK4_INTERVALLISTTOOLS_SCATTER } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: ['args':'--SCATTER_COUNT 6'] )
include { GATK4_INTERVALLISTTOOLS AS GATK4_INTERVALLISTTOOLS_SUBDIVISION } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: ['args':'--SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW'] )
include { GATK4_INTERVALLISTTOOLS AS GATK4_INTERVALLISTTOOLS_UNIQUE } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: ['args':'--UNIQUE true'] )
include { GATK4_INTERVALLISTTOOLS AS GATK4_INTERVALLISTTOOLS_SORT } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: ['args':'--SORT true'] )


workflow test_gatk4_intervallisttools_scatter {

    input = [ [ id:'test' ], // meta map
              file('/home/praveen/Desktop/tools/downloads/human_gtf/human.exons.interval_list', checkIfExists: true) ]

    GATK4_INTERVALLISTTOOLS_SCATTER ( input )
}


workflow test_gatk4_intervallisttools_subdivision {

    input = [ [ id:'test' ], // meta map
              file('/home/praveen/Desktop/tools/downloads/human_gtf/human.exons.interval_list', checkIfExists: true) ]

    GATK4_INTERVALLISTTOOLS_SUBDIVISION ( input )
}


workflow test_gatk4_intervallisttools_unique {

    input = [ [ id:'test' ], // meta map
              file('/home/praveen/Desktop/tools/downloads/human_gtf/human.exons.interval_list', checkIfExists: true) ]

    GATK4_INTERVALLISTTOOLS_UNIQUE ( input )
}


workflow test_gatk4_intervallisttools_sort {

    input = [ [ id:'test' ], // meta map
              file('/home/praveen/Desktop/tools/downloads/human_gtf/human.exons.interval_list', checkIfExists: true) ]

    GATK4_INTERVALLISTTOOLS_SORT ( input )
}
