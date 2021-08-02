#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RMARKDOWN } from '../../../modules/rmarkdown/main.nf' addParams(
     options: [parametrize: false]
)
include { RMARKDOWN_PARAMETRIZE } from '../../../modules/rmarkdown/main.nf' addParams(
    options: [:]
)

workflow test_rmarkdown {

    input = [ [ id:'test_rmd' ], // meta map
              file(params.test_data['generic']['notebooks']['rmarkdown'], checkIfExists: true) ]

    RMARKDOWN ( input, [:], [:])

}

workflow test_rmarkdown_parametrize {

    input = [ [ id:'test_rmd' ], // meta map
              file(params.test_data['generic']['notebooks']['rmarkdown'], checkIfExists: true) ]

    RMARKDOWN_PARAMETRIZE(
        input,
        [input_filename: "hello.txt", n: 12],
        file(params.test_data['generic']['txt']['hello'])
    )

}

