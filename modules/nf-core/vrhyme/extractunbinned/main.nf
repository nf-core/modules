// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VRHYME_EXTRACTUNBINNED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e2b5a580b1eb2e3b5d4c50ba78df05d4e49e0bb1dfc988e719306fa9efb21d65/data'
:         'community.wave.seqera.io/library/vrhyme_python_setuptools:2e133279ed381753' }"

    input:
    tuple val(meta), path(membership)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fasta") , emit: unbinned_sequences
    tuple val("${task.process}"), val('vrhyme'), eval("vRhyme --version 2>&1 | sed 's/^.*vRhyme v//; s/Using.*\$//'"), emit: versions_vrhyme, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_unbinned_sequences"
    def fasta_input = fasta.toString().replaceAll(/\.gz$/, '')
    def gunzip      = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${fasta_input}" : ""
    def cleanup     = fasta.getExtension() == "gz" ? "rm ${fasta_input}" : ""
    if( "$fasta_input" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    ${gunzip}

    extract_unbinned_sequences.py \\
        -i ${membership} \\
        -f ${fasta_input} \\
        -o ${prefix}.fasta \\
        ${args}

    ${cleanup}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_unbinned_sequences"
    def fasta_input = fasta.toString() - ~/\.gz$/
    if( "$fasta_input" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta
    """
}
