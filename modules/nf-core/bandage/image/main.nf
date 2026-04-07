process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3eabbd074e3bc45e2643783450330cae3afc6697fefc635755ab964dc43665a1/data' :
        'community.wave.seqera.io/library/bandage:0.9.0--4f0567049a14ea6d' }"

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

    Bandage image ${gfa_input} ${prefix}.png $args
    Bandage image ${gfa_input} ${prefix}.svg $args

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}.svg
    """
}
