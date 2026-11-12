// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process OPENMSTHIRDPARTY_PERCOLATORADAPTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/6475e53b3776cd3b71b8ff525f01da9c45dbe9401a00e5630d91ab976f20e8aa/data'
        : 'community.wave.seqera.io/library/openms-thirdparty:3.6.0--566f6f290d39203b'}"

    input:
    tuple val(meta), path(id_file)

    output:
    tuple val(meta), path("${prefix}.idXML"), emit: idxml
    tuple val(meta), path("*_feature_weights.tsv"), emit: feature_weights, optional: true
    tuple val("${task.process}"), val('PercolatorAdapter'), eval("PercolatorAdapter --help 2>&1 | grep -oP '^Version: \\K[0-9.]+'"), emit: versions_percolatoradapter, topic: versions
    tuple val("${task.process}"), val('percolator'), eval("percolator -h 2>&1 | grep -oP 'version \\K[0-9.]+'"), emit: versions_percolator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_percolator"
    if ("${id_file}" == "${prefix}.idXML") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    PercolatorAdapter \\
        -in ${id_file} \\
        -out ${prefix}.idXML \\
        -threads ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_percolator"
    if ("${id_file}" == "${prefix}.idXML") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    echo ${args}

    touch ${prefix}.idXML
    touch ${prefix}_feature_weights.tsv
    """
}
