#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'

workflow test_ampcombi {
    file_list = channel.fromPath(
       [ 
            'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/ampir/sample_2/sample_2.ampir.tsv',
            'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/amplify/sample_2/sample_2.amplify.tsv',
        ] , checkIfExists: true ).collect()
    
    samplesheet = file("https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/samplesheet.csv", checkIfExists: true)

    faa_folder = file("https://github.com/nf-core/test-datasets/tree/modules/data/delete_me/ampcombi/test_faa/", checkIfExists: true)
    
    outdir = "ampcombi_results"

    AMPCOMBI ( file_list, samplesheet, faa_folder, outdir )
}

//['https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/ampir/sample_2/sample_2.ampir.tsv',
//'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/amplify/sample_2/sample_2.amplify.tsv']
//file("https://raw.githubusercontent.com/nf-core/test-datasets/funcscan/samplesheet.csv", checkIfExists: true)                                                                                                                            
//"\$SAMPLE"