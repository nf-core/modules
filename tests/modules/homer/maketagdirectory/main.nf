#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKETAGDIRECTORY } from '../../../../software/homer/maketagdirectory/main.nf' addParams( options: [:] )

workflow test_homer_annotatepeaks {
    HOMER_MAKETAGDIRECTORY (  )
}

