// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process OPENMS_IDMASSACCURACY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/233b3693133c970e693ef74e8a599dd8984d987e4d3e298f43bf520dd184f754/data'
:         'community.wave.seqera.io/library/openms:3.4.1--6825b96bd18a1e9e' }"

    input:
    tuple val(meta), path(mzmls), path(idxmls)

    output:
    tuple val(meta), path("*frag_mass_err.tsv") , emit: frag_err
    tuple val(meta), path("*prec_mass_err.tsv") , emit: prec_err, optional: true
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

    Reason:
    This module is no longer fit for purpose because not part of openms 3.5.0 version anymore
    """
    assert false: deprecation_message

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    IDMassAccuracy \\
        -in $mzmls \\
        -id_in $idxmls \\
        -out_fragment ${prefix}_frag_mass_err.tsv \\
        -threads $task.cpus \\
        $args

    """

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

    Reason:
    This module is no longer fit for purpose because not part of openms 3.5.0 version anymore
    """
    assert false: deprecation_message
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_frag_mass_err.tsv
    touch ${prefix}_prec_mass_err.tsv

    """
}
