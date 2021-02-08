#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_FINDPEAKS } from '../../../../software/homer/findpeaks/main.nf' addParams( options: [:] )

workflow test_homer_annotatepeaks {
    HOMER_FINDPEAKS (  )
}

