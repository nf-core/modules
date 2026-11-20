// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59db1f1331e9bd80560221e6283dea4e9f14de41e2f028c13b71f0b9e42b8afb/data' :
        'community.wave.seqera.io/library/bismark:3.3.0--f5ac4e7eb5a16146' }"

    input:
    tuple val(meta), path(coverage_file)
    tuple val(meta2), path(fasta, stageAs: 'tmp/*') // This change mounts as directory containing the FASTA file to prevent nested symlinks
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.cov.gz")                      , emit: coverage,  optional: true
    tuple val(meta), path("*report.txt.gz")                , emit: report
    tuple val(meta), path("*cytosine_context_summary.txt") , emit: summary
    tuple val("${task.process}"), val('bismark'), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), emit: versions_bismark, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage2cytosine \\
        ${coverage_file} \\
        --genome ${index} \\
        --output ${prefix} \\
        --gzip \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.cov.gz
    echo "" | gzip > ${prefix}.report.txt.gz
    touch ${prefix}.cytosine_context_summary.txt
    """
}
