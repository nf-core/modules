#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABULI_BUILD } from '../../../../../modules/nf-core/metabuli/build/main.nf'

process CREATE_TAXONOMY_FOLDER{
    input:
    path(dmpfiles)

    output:
    path(db)

    script:
    """
    mkdir -p db/taxonomy
    mv *.dmp db/taxonomy
    touch db/taxonomy/merged.dmp
    """

}

workflow test_metabuli_build {
    
    genome = file("${params.test_data_base}/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    
    dmp_files = [
        file("${params.test_data_base}/data/delete_me/metabuli/taxonomy/names.dmp"),
        file("${params.test_data_base}/data/delete_me/metabuli/taxonomy/nodes.dmp")
    ]
    acc2taxid = file("${params.test_data_base}/data/delete_me/metabuli/acc2taxid")

    tax = CREATE_TAXONOMY_FOLDER(dmp_files) 
    METABULI_BUILD ( genome, acc2taxid, tax )
}
