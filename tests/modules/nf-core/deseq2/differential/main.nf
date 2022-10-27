#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { DESEQ2_DIFFERENTIAL } from "$moduleDir/modules/nf-core/deseq2/differential/main.nf"

process tsv_to_csv {

    input:
    path expression_matrix

    output:
    path 'test.csv'

    script:
    """
    sed 's/\t/,/g' ${expression_matrix} > test.csv.tmp
    mv test.csv.tmp test.csv
    """

}

workflow test_deseq2_differential {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, expression_sample_sheet, expression_matrix_file)
        }
        .set{
            input
        }

    DESEQ2_DIFFERENTIAL (
        input
    )
}

// This second test checks that a CSV format matrix works the same as a TSV one

workflow test_deseq2_differential_csv {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)


    ch_samples_and_matrix  = tsv_to_csv(expression_matrix_file)
        .map{
            tuple(expression_sample_sheet, it)
        }

    ch_contrasts = Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )

    DESEQ2_DIFFERENTIAL (
        ch_contrasts.combine(ch_samples_and_matrix)
    )
}
