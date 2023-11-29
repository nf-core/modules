#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_TABULARTOGSEAGCT } from '../../../../../modules/nf-core/custom/tabulartogseagct/main.nf'

infile = params.test_data['mus_musculus']['genome']['rnaseq_matrix']

workflow test_custom_tabulartogseagct {

    input = Channel.fromPath(infile)
        .splitCsv(sep: "\t", header: false)
        .map{
            lst = new ArrayList(it);
            lst.remove(1);
            lst.join('\t')
        }
        .collectFile(name: 'test.tsv', newLine: true, sort: false)
        .map{
            [ [ id:'test' ], it]
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}

workflow test_custom_tabulartogseagct_csv {
    
    input = Channel.fromPath(infile)
        .splitCsv(sep: "\t", header: false)
        .map{
            lst = new ArrayList(it);
            lst.remove(1);
            lst.join(',')
        }
        .view()
        .collectFile(name: 'test.csv', newLine: true, sort: false)
        .map{
            [ [ id:'test' ], it]
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}

workflow test_custom_tabulartogseagct_csv_override {
    
    input = Channel.fromPath(infile)
        .splitCsv(sep: "\t", header: false)
        .map{
            lst = new ArrayList(it);
            lst.remove(1);
            lst.join(',')
        }
        .view()
        .collectFile(name: 'test.tsv', newLine: true, sort: false)
        .map{
            [ [ id:'test' ], it]
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}

workflow test_custom_tabulartogseagct_csv_override_pipe {
    
    input = Channel.fromPath(infile)
        .splitCsv(sep: "\t", header: false)
        .map{
            lst = new ArrayList(it);
            lst.remove(1);
            lst.join('|')
        }
        .view()
        .collectFile(name: 'test.tsv', newLine: true, sort: false)
        .map{
            [ [ id:'test' ], it]
        }

    CUSTOM_TABULARTOGSEAGCT ( input )
}
