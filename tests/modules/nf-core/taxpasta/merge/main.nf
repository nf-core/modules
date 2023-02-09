#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TAXPASTA_MERGE } from '../../../../../modules/nf-core/taxpasta/merge/main.nf'


workflow test_taxpasta_merge {


test1_kraken = Channel.fromPath(params.test_data['sarscov2']['metagenome']['kraken_report']).collectFile( name: 'test1_kraken')
test2_kraken = Channel.fromPath(params.test_data['sarscov2']['metagenome']['kraken_report']).collectFile( name: 'test2_kraken')

input = test1_kraken.mix ( test2_kraken ).collect()
                    .map { files ->
                               def meta = [:]
                               meta['id'] = 'kraken2'
                               meta['profiler'] = 'kraken2'
                          [meta,files]
                    }


	TAXPASTA_MERGE ( input, [] )
}


