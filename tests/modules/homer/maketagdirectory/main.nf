#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../software/homer/maketagdirectory/main.nf' addParams( options: [:] )

    HOMER_MAKETAGDIRECTORY (  )
workflow test_homer_maketagdirectory {
}

