process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3eabbd074e3bc45e2643783450330cae3afc6697fefc635755ab964dc43665a1/data' :
        'community.wave.seqera.io/library/bandage:0.9.0--4f0567049a14ea6d' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Bandage image $gfa ${prefix}.png $args
    Bandage image $gfa ${prefix}.svg $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(echo \$(export QT_QPA_PLATFORM=offscreen; Bandage --version 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
