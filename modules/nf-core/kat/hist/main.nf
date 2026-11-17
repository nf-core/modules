// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KAT_HIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b722aba17c3463207a656322aab4c4a443f904fef5b305809eed489a2412489/data'
        : 'community.wave.seqera.io/library/kat:2.4.2--75e8f5922e348f7c'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hist")                   , emit: hist
    tuple val(meta), path("*.hist.dist_analysis.json"), emit: json
    tuple val(meta), path("*.png")                    , emit: png           , optional: true
    tuple val(meta), path("*.ps")                     , emit: ps            , optional: true
    tuple val(meta), path("*.pdf")                    , emit: pdf           , optional: true
    tuple val(meta), path("*-hash.jf*")               , emit: jellyfish_hash, optional: true
    tuple val("${task.process}"), val('kat'), eval("kat hist --version | sed 's/kat //'"), emit: versions_kat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/merqury/merqury

Reason:
This module no longer works in conda due to glibc incompatibilities with plotting libraries
This module is no longer maintained by the authors
"""
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    kat hist \\
        --threads $task.cpus \\
        --output_prefix ${prefix}.hist \\
        $args \\
        $reads

    ls -l

    """

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/merqury/merqury

Reason:
This module no longer works in conda due to glibc incompatibilities with plotting libraries
This module is no longer maintained by the authors
"""
    def prefix    = task.ext.prefix ?: "${meta.id}"
    assert false: deprecation_message
    """
    touch ${prefix}.hist

    """
}
