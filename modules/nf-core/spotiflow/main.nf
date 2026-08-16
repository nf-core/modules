process SPOTIFLOW {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f7e014ce832214c25d8faa1b4660db125a79c8c3b075c93059771aa3b83d1c3/data':
        'community.wave.seqera.io/library/spotiflow:0.5.7--c7eb617591164164' }"

    input:
    tuple val(meta), path(image_2d)

    output:
    tuple val(meta), path("*.csv"), emit: spots
    tuple val("${task.process}"), val('spotiflow'), eval("python -m pip show spotiflow | grep 'Version' | sed -e 's/Version: //g'"), emit: versions_spotiflow, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    """
    spotiflow-predict \\
        ${image_2d} \\
        --out-dir . \\
        $args
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
