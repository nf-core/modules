#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { JUPYTERNOTEBOOK } from '../../../modules/jupyternotebook/main.nf'
include { JUPYTERNOTEBOOK as JUPYTERNOTEBOOK_PARAMETRIZE } from '../../../modules/jupyternotebook/main.nf'
include { JUPYTERNOTEBOOK as JUPYTERNOTEBOOK_PARAMETRIZE_IPYNB } from '../../../modules/jupyternotebook/main.nf'

workflow test_jupyternotebook {

    input = [ [ id:'test_jupyter' ], // meta map
              file(params.test_data['generic']['notebooks']['ipython_md'], checkIfExists: true) ]

    JUPYTERNOTEBOOK ( input, [:], [])

}

workflow test_jupyternotebook_parametrize {

    input = [ [ id:'test_jupyter' ], // meta map
              file(params.test_data['generic']['notebooks']['ipython_md'], checkIfExists: true) ]

    JUPYTERNOTEBOOK_PARAMETRIZE(
        input,
        [input_filename: "hello.txt", n_iter: 12],
        file(params.test_data['generic']['txt']['hello'], checkIfExists: true)
    )

}

workflow test_jupyternotebook_parametrize_ipynb {

    input = [ [ id:'test_jupyter' ], // meta map
              file(params.test_data['generic']['notebooks']['ipython_ipynb'], checkIfExists: true) ]

    JUPYTERNOTEBOOK_PARAMETRIZE_IPYNB(
        input,
        [input_filename: "hello.txt", n_iter: 12],
        file(params.test_data['generic']['txt']['hello'], checkIfExists: true)
    )

}

