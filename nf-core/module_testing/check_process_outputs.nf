#!/usr/bin/env nextflow
nextflow.preview.dsl = 2

cheers = Channel.from 'Bonjour', 'Ciao', 'Hello', 'Hola'

process check_output {
    input:
    val x from cheers

    script:
    """
    echo '$x world!'
    """
}
