process CLAME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clame:1.0--he1b5a44_1':
        'biocontainers/clame:1.0--he1b5a44_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta")  , emit: fasta, optional: true
    tuple val(meta), path("*.binning"), emit: bins
    tuple val(meta), path("*.fm9")    , emit: fm
    tuple val(meta), path("*.index")  , emit: index
    tuple val(meta), path("*.links")  , emit: links
    tuple val(meta), path("*.result") , emit: result
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clame \\
        $args \\
        -nt $task.cpus \\
        -multiFasta ${fasta} \\
        -output ${prefix} || test -f ${prefix}.binning
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clame: \$(echo \$(clame -h | sed -n '2p' | cut -d ' ' -f 2 ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.binning
    touch ${prefix}.fm9
    touch ${prefix}.index
    touch ${prefix}.links
    touch ${prefix}.result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clame: \$(echo \$(clame -h | sed -n '2p' | cut -d ' ' -f 2 ))
    END_VERSIONS
    """
}
