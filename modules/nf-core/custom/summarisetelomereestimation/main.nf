// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_SUMMARISETELOMEREESTIMATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/93f76f0461c25ad04c5f645d7eadbc8ec12c7b2100b0ffef05dc6509e03dac3e/data'
        : 'community.wave.seqera.io/library/pandas_python:a1229de030b588d2'}"

    input:
    tuple val(meta), path(length_tsv), path(content_tsv), val(length_tool)

    output:
    tuple val(meta), path("*_telomere_summary.tsv"), emit: summary
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'summarisetelomereestimation.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_telomere_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
