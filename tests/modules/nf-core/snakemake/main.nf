#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNAKEMAKE } from '../../../../modules/nf-core/snakemake/main.nf'

workflow test_snakemake {

    input = [
        [ id: 'input'],
        []
    ]

// This generates a Snakefile for use with Snakemake.
// In real use, you should probably access this via a
// a file stored within the repository.
    Channel.of('''
rule all:
    input: "hello.txt"

rule hello_world:
    output: "hello.txt"
    shell: "echo Hello World > hello.txt"

'''
        )
        .collectFile(name: 'Snakefile')
        .map { file ->
            [
                [id: 'Snakefile'],
                file
            ]}
        .set{ snakefile }
    snakefile.view()

    SNAKEMAKE ( input, snakefile)
}
