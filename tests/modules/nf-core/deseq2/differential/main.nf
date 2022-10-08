#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.deseq_vs_method = 'rlog'
params.deseq_random_seed = 1234
params.deseq_round_results = 'TRUE'

include { DESEQ2_DIFFERENTIAL } from '../../../../../modules/nf-core/deseq2/differential/main.nf'

process spoof_samplesheet {

    input:
    path expression_matrix
    val contrast
    
    output:
        path 'samplesheet.csv' 

    script:
    """
    echo "experiment_accession,treatment" > tmp.csv
    head -n 1 $expression_matrix | tr "\t" "\n" | tail -n 6 | head -n 3 | while read -r l; do
        echo \$l,$contrast.reference >> tmp.csv
    done 
    head -n 1 $expression_matrix | tr "\t" "\n" | tail -n 3 | while read -r l; do
        echo \$l,$contrast.target >> tmp.csv
    done 
    
    mv tmp.csv samplesheet.csv
    """

}

workflow test_deseq2_differential {

    ch_contrast_definitions = Channel.from([ "variable":"treatment", "reference":"saline", "target":"drug", "blocking":"" ])

    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 

    DESEQ2_DIFFERENTIAL (
        spoof_samplesheet(
            expression_matrix_file,
            ch_contrast_definitions
        ),
        expression_matrix_file,
        ch_contrast_definitions     
    )
}
