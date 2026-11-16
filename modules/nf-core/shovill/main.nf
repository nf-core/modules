// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SHOVILL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/97e3f0d4872793310b8bada5c56aa3905bb2ac430b91a7521a84a06b72e1768f/data'
        : 'community.wave.seqera.io/library/shovill:1.4.2--eab27bd761ae82ef'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("contigs.fa")                         , emit: contigs
    tuple val(meta), path("shovill.corrections")                , emit: corrections
    tuple val(meta), path("shovill.log")                        , emit: log
    tuple val(meta), path("{skesa,spades,megahit,velvet}.fasta"), emit: raw_contigs
    tuple val(meta), path("{skesa,spades,megahit,velvet}.{fastg,gfa,LastGraph}"), optional:true, emit: gfa
    tuple val("${task.process}"), val('shovill'), eval('shovill --version 2>&1 | sed "s/^.*shovill //"'), emit: versions_shovill, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory.toGiga()
    """
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $args \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force
    """

    stub:
    """
    touch contigs.fa
    touch shovill.corrections
    touch shovill.log
    touch spades.fasta
    """
}
