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
ch_empty_spikes = [[],[]]

workflow test_deseq2_differential {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
    )
}

// Not including blocking column in contrasts file shouldn't be lethal

workflow test_deseq2_differential_noblocking {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            it.remove('blocking')
            tuple(it, it.variable, it.reference, it.target)
        }
        .unique()
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
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
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    // Make our fake spikes and pretend they're ERCC controls

    spoof_spikes(expression_matrix_file)
        .map{
            tuple(['id':'ERCC'], it)
        }.set{
            ch_spikes
        }

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
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
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    // Make our fake spikes and pretend they're ERCC controls

    spoof_spikes(expression_matrix_file)
        .map{
            tuple(['id':'ERCC'], it)
        }.set{
            ch_spikes
        }

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_spikes
    )
}

// This second test checks that a CSV format matrix works the same as a TSV one

workflow test_deseq2_differential_csv {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    ch_matrix  = tsv_to_csv(expression_matrix_file)
        .map{
            tuple([id: 'test'], expression_sample_sheet, it)
        }

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    DESEQ2_DIFFERENTIAL(
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
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
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
    )
}

// Test with contrast sample subsetting enabled

workflow test_deseq2_differential_subset_to_contrast {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
    )
}

// Test with sample exclusion enabled

workflow test_deseq2_differential_exclude_samples {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true)
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true)
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true)

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(it, it.variable, it.reference, it.target)
        }
        .set{
            ch_contrasts
        }

    ch_matrix = [[id: 'test'], expression_sample_sheet, expression_matrix_file]

    DESEQ2_DIFFERENTIAL (
        ch_contrasts,
        ch_matrix,
        ch_empty_spikes
    )
}
