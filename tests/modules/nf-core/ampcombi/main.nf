#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPCOMBI } from '../../../../modules/nf-core/ampcombi/main.nf'

workflow test_ampcombi {
    file_list = "['https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/ampir/sample_1/sample_1.ampir.tsv','https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/amplify/sample_1/sample_1_amplify.tsv']['https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/ampir/sample_2/sample_2.ampir.tsv','https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/ampcombi/test_files/amplify/sample_2/sample_2.amplify.tsv']"

    faa_folder = "https://github.com/nf-core/test-datasets/tree/modules/data/delete_me/ampcombi/test_faa"
    
    outdir = "ampcombi_results"

    AMPCOMBI ( file_list, faa_folder, outdir )
}