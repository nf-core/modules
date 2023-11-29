#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANNOTSV_INSTALLANNOTATIONS } from '../../../../../modules/nf-core/annotsv/installannotations/main.nf'

workflow test_annotsv_installannotations {    
    ANNOTSV_INSTALLANNOTATIONS()
}
