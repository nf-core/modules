#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAGECK_COUNT } from '../../../../modules/mageck/count/main.nf'

workflow test_mageck_count_fastq {
    def input = []
    input = [ [ id:'test,test2', single_end:true] , // meta map
        [file('https://raw.githubusercontent.com/nf-core/test-datasets/crisprseq/testdata/ERR376998.small.fastq.gz', checkIfExists: true),
        file('https://raw.githubusercontent.com/nf-core/test-datasets/crisprseq/testdata/ERR376999.small.fastq.gz', checkIfExists: true)]
    ]
    library = file('https://raw.githubusercontent.com/nf-core/test-datasets/crisprseq/testdata/yusa_library.csv', checkIfExists: true)

    MAGECK_COUNT ( input, library )
}


workflow test_mageck_count_counts {
    def input = []
    input =  [[ id:'test', single_end:true] , // meta map
        [file('https://raw.githubusercontent.com/nf-core/test-datasets/crisprseq/testdata/count_table.csv', checkIfExists: true)]]
    library = file('https://raw.githubusercontent.com/nf-core/test-datasets/crisprseq/testdata/yusa_library.csv', checkIfExists: true)


    MAGECK_COUNT ( input, library )
}
