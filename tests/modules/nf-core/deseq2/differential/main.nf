#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DESEQ2_DIFFERENTIAL } from '../../../../../modules/nf-core/deseq2/differential/main.nf'

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

// Take the first 50 genes and pretend they're spikes for testing control gene
// functionality

process spoof_spikes {

    input:
    path expression_matrix

    output:
    path 'spikes.txt'

    script:
    """
    head -n 50 $expression_matrix | \
        tail -n +2 | \
        awk '{print \$1}' > spikes.txt.tmp
    mv spikes.txt.tmp spikes.txt
    """
}
empty_spikes = [[],[]]

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
        input,
        empty_spikes
    )
}

// Try with spikes as control genes

workflow test_deseq2_differential_spikes {

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

    // Make our fake spikes and pretend they're ERCC controls

    spoof_spikes(expression_matrix_file)
        .map{
            tuple(['id':'ERCC'], it)
        }.set{
            ch_spikes
        }

    DESEQ2_DIFFERENTIAL (
        input,
        ch_spikes
    )
}

// Try with spikes as control genes, but stripping rather than using

workflow test_deseq2_differential_strip_spikes {

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

    // Make our fake spikes and pretend they're ERCC controls

    spoof_spikes(expression_matrix_file)
        .map{
            tuple(['id':'ERCC'], it)
        }.set{
            ch_spikes
        }

    DESEQ2_DIFFERENTIAL (
        input,
        ch_spikes
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
        ch_contrasts.combine(ch_samples_and_matrix),
        empty_spikes
    )
}

// Test vst with modified nsub

workflow test_deseq2_differential_vst_nsub {

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
        input,
        empty_spikes
    )
}
