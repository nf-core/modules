process GECCO_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gecco:0.10.1--pyhdfd78af_0':
        'biocontainers/gecco:0.10.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(clusters), path(gbk)
    val(mode)
    val(format)

    output:
    tuple val(meta), path("${prefix}/*.gff")        , emit: gff     , optional: true
    tuple val(meta), path("${prefix}/*.region*.gbk"), emit: bigslice, optional: true
    tuple val(meta), path("${prefix}/*.faa")        , emit: faa     , optional: true
    tuple val(meta), path("${prefix}/*.fna")        , emit: fna     , optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gecco \\
        convert \\
        $args \\
        --jobs $task.cpus \\
        $mode \\
        --input-dir ./ \\
        --format ${format} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: \$(echo \$(gecco --version) | cut -f 2 -d ' ' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir ${prefix}
    touch ${prefix}/${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gecco: \$(echo \$(gecco --version) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
