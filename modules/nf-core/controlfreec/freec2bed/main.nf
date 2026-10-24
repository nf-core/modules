// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CONTROLFREEC_FREEC2BED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/86e7bfef8a43af72880501e04fb817268228715f4854102ce4a90fd3ad2cd66a/data'
:         'community.wave.seqera.io/library/control-freec:11.6b--6482bd27b629122a' }"

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('controlfreec'), val("11.6b"), emit: versions_controlfreec, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freec2bed.pl -f ${ratio} ${args} > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
