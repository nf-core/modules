process EMBOSS_GETORF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5':
        'biocontainers/emboss:6.6.0--h86d058a_5' }"

    input:
    tuple val(meta), path(sequence)
    val out_ext

    output:
    tuple val(meta), path("*.${out_ext}"), emit: orf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def osformat2 = args.contains('-osformat2') ? '' : "-osformat2 ${out_ext}"
    def table = params.getorf_table ?: 0
    def minsize = params.getorf_minsize ?: 0
    def maxsize = params.getorf_maxsize ?: 1000000
    def getorf_find = params.getorf_find ?: 0

    """
    getorf \\
    -table ${table} \\
    -minsize ${minsize} \\
    -maxsize ${maxsize} \\
    -find ${getorf_find} \\
    -sequence ${sequence} \\
    ${osformat2} \\
    -outseq ${prefix}.${out_ext} \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(echo \$(getorf -version 2>&1) | sed 's/EMBOSS://')
    END_VERSIONS
    """
}
