#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FREYJA_UPDATE } from '../../../../../modules/nf-core/freyja/update/main.nf'

workflow test_freyja_update {
    db_name= "freyja_db"
    FREYJA_UPDATE (db_name)
}
