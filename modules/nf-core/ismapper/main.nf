// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ISMAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8db3e77ec9865b4d89d154d39f4d918a9eb1d91b446f92d21d836634b98d682/data'
        : 'community.wave.seqera.io/library/ismapper_python_setuptools_samtools:e034790cb6b187c4'}"

    input:
    tuple val(meta), path(reads), path(reference), path(query)

    output:
    tuple val(meta), path("results/*"), emit: results
    tuple val("${task.process}"), val('ismapper'), eval("ismap --version 2>&1 | sed -n 's/^.*ismap //p'"), emit: versions_ismapper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ismap \\
        ${args} \\
        --t ${task.cpus} \\
        --output_dir results \\
        --queries ${query} \\
        --reference ${reference} \\
        --reads ${reads}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/${prefix}

    touch results/${prefix}/${prefix}_left_final.fastq
    """
}
