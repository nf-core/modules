process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24e7692f7cb34d64c5e370dc9dac2e2b6f687d37b7140f12a6251d25dd8b6fa6/data'
        : 'community.wave.seqera.io/library/bandage:0.9.0--5aa7953644fadc5e'}"

    input:
    tuple val(meta), path(gfa, arity: '1')

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    tuple val("${task.process}"), val('bandage'), eval('export QT_QPA_PLATFORM=offscreen; Bandage --version 2>&1| grep Version | sed "s/^Version: //"'), emit: versions_bandage, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def gfa_input = gfa.toString().endsWith('.gz') ? gfa.toString() - ~/\.gz$/ : gfa.toString()
    def decompress = gfa.toString().endsWith('.gz') ? "zcat ${gfa} > ${gfa_input}" : ""
    def cleanup = gfa.toString().endsWith('.gz') ? "rm ${gfa_input}" : ""
    """
    ${decompress}

    Bandage image ${gfa_input} ${prefix}.png ${args}
    Bandage image ${gfa_input} ${prefix}.svg ${args}

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}.svg
    """
}
