#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_CONFIGUREHOMER } from '../../../../software/homer/configurehomer/main.nf' addParams( options: [:] )

workflow test_homer_configurehomer {
      HOMER_CONFIGUREHOMER ( )
}

