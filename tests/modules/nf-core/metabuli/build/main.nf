#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABULI_BUILD } from '../../../../../modules/nf-core/metabuli/build/main.nf'

process BUILD_ACC2TAXID {
    input:
      path(genomes)
    output:
      path(acc2taxid)

    script:
    """
    echo -e "accession\taccession.version\ttaxid\tgi" > acc2taxid
    accessionv=\$(cat ${genomes} | head -n1 | 
                  cut -d " " -f1 | sed 's/>//')
    echo -e "\${accessionv%.*}\t\$accessionv\t2697049\t111" >> acc2taxid
    """
}

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
        file("${params.test_data_base}/delete_me/metabuli/names.dmp"),
        file("${params.test_data_base}/delete_me/metabuli/nodes.dmp")
    ]
    acc2taxid = BUILD_ACC2TAXID(genome)
    tax = CREATE_TAXONOMY_FOLDER(dmp_files) 
    METABULI_BUILD ( genome, acc2taxid, tax )
}
