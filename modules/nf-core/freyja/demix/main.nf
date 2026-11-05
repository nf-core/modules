// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FREYJA_DEMIX {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f0b023f9c9e25f338d72ed6710e4d2a54a4c3a9a3ac34f17e08683894dc40830/data'
:         'community.wave.seqera.io/library/freyja:2.0.6--be3cbd6f14f9d542' }"

    input:
    tuple val(meta), path(variants), path(depths)
    path barcodes
    path lineages_meta
    path lineages_topology

    output:
    tuple val(meta), path("*.tsv"), emit: demix
    tuple val("${task.process}"), val('freyja'), eval("freyja --version | sed 's/.* //'"), topic: versions, emit: versions_freyja

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meta_cmd = lineages_meta ? "--meta $lineages_meta" : ''
    def lineage_cmd = lineages_topology ? "--lineageyml $lineages_topology" : ''
    """
    freyja \\
        demix \\
        $args \\
        --output ${prefix}.tsv \\
        --barcodes $barcodes \\
        $lineage_cmd \\
        $meta_cmd \\
        $variants \\
        $depths
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    """

}
