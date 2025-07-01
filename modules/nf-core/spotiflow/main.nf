process SPOTIFLOW {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/spotiflow:77826b49fcab1875':
        'community.wave.seqera.io/library/spotiflow:b300f293b176b30e' }"

    input:
    tuple val(meta), path(image_2d)

    output:
    tuple val(meta), path("*.csv"), emit: spots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    """
    spotiflow-predict \\
        ${image_2d} \\
        --out-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spotiflow: \$( python -m pip show --version spotiflow | grep "Version" | sed -e "s/Version: //g" )
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spotiflow: \$( python -m pip show --version spotiflow | grep "Version" | sed -e "s/Version: //g" )
    END_VERSIONS
    """
}
