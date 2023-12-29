process MUDSKIPPER_INDEX {
    tag '$gtf'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mudskipper:0.1.0--h9f5acd7_1':
        'biocontainers/mudskipper:0.1.0--h9f5acd7_1' }"

    input:
    path gtf

    output:
    path "index/"      , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export RUST_BACKTRACE=full
    mudskipper \\
        index \\
        --gtf $gtf \\
        --dir-index index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mudskipper: \$(echo \$(mudskipper -V 2>&1) | sed 's/^.*mudskipper //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir index/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mudskipper: \$(echo \$(mudskipper -V 2>&1) | sed 's/^.*mudskipper //' )
    END_VERSIONS
    """
}
