process SEQKIT_SANA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.1--he881be0_0':
        'biocontainers/seqkit:2.10.1--he881be0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}${extension}"), emit: reads
    tuple val(meta), path("${prefix}.log")        , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = reads.getName() - reads.getSimpleName()
    """
    seqkit sana \\
        $args \\
        ${reads} \\
        -o ${prefix}${extension} > ${prefix}.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = reads.getName() - reads.getSimpleName()
    """
    echo $args

    # Create empty output files conditionally based on extension
    if [[ "${extension}" == *.gz ]]; then
        echo -n | gzip -c > ${prefix}${extension}
    else
        touch ${prefix}${extension}
    fi
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
