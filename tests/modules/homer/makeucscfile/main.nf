#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_MAKEUCSCFILE } from '../../../../software/homer/makeucscfile/main.nf' addParams( options: [:] )

workflow test_homer_annotatepeaks {
    HOMER_MAKEUCSCFILE (  )
}

