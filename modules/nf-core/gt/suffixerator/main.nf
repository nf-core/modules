process GT_SUFFIXERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("$prefix"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir \\
        "$prefix"

    gt \\
        suffixerator \\
        $args \\
        -db $fasta \\
        -indexname "$prefix/suffixerator"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir \\
        "$prefix"

    touch "$prefix/suffixerator.esq"
    touch "$prefix/suffixerator.prj"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
