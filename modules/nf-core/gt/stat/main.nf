process GT_STAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("${prefix}.yml")  , emit: stats
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    gt \\
        stat \\
        $args \\
        $gff3 \\
        > ${prefix}.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
