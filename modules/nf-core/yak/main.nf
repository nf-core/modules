process YAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::yak=0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yak:0.1--h7132678_2':
        'biocontainers/yak:0.1--h7132678_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.yak"), emit: yak
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    yak \\
        count \\
        $args \\
        -t${task.cpus} \\
        -o ${prefix}.yak <(zcat $reads)


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak : \$(echo \$(yak --version 2>&1) | sed 's/^.*yak //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yak

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak : \$(echo \$(yak --version 2>&1) | sed 's/^.*yak //' )
    END_VERSIONS
    """
}
