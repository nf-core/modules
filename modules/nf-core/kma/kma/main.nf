process KMA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.15--h577a1d6_1' :
        'biocontainers/kma:1.4.15--h577a1d6_1' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    val interleaved
    val mat_format

    output:
    tuple val(meta), path("*.res")      , optional: true, emit: res
    tuple val(meta), path("*.fsa")      , optional: true, emit: fsa
    tuple val(meta), path("*.aln")      , optional: true, emit: aln
    tuple val(meta), path("*.frag.gz")  , optional: true, emit: frag
    tuple val(meta), path("*.mat.gz")   , optional: true, emit: mat   // if mat_format == true
    tuple val(meta), path("*.spa")      , optional: true, emit: spa   // if ext.args contains '-Sparse' (only output in this case)
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}.${meta2.id}"
    def input_style     = interleaved ? "-int ${reads}" : "-ipe ${reads}"
    def create_mat      = mat_format ? "-matrix" : ''
    """
    INDEX=`find -L ./ -name "*.name" | sed 's/\\.name\$//'`

    kma \\
        ${input_style} \\
        -o ${prefix} \\
        -t_db \$INDEX \\
        ${create_mat} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-\$//')
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta2.id}"
    def create_spa  = task.ext.args.contains('-Sparse')

    if ( create_spa )
        """
        touch ${prefix}.spa
        """
    else
        """
        touch ${prefix}.res \\
        touch ${prefix}.fsa \\
        touch ${prefix}.aln \\
        echo "" | gzip > ${prefix}.frag.gz
        """

    if ( mat_format )
        """
        echo "" | gzip > ${prefix}.mat.gz
        """

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v 2>&1) | sed 's/^KMA-\$//')
    END_VERSIONS
    """
}
