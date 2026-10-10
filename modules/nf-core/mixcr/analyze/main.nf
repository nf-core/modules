process MIXCR_ANALYZE {
    tag "$meta.id"
    label 'process_high'

    secret 'MI_LICENSE'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7dccca2db544d708a462b2ba49117b41e32141d8eb46f4dbc8730e210d2cdbae/data':
        'community.wave.seqera.io/library/mixcr:4.7.0--bb57944ca92aeb74' }"

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        '-B \$HOME' :
        '' }"
    
    input:
    tuple val(meta), path(reads)
    val preset
    val species

    output:
    tuple val(meta), path("*clones*.tsv")    , emit: clones
    tuple val(meta), path("*.txt")           , emit: reports, optional: true
    tuple val(meta), path("*.contigs.clns")  , emit: clns, optional: true
    tuple val(meta), path("*.extended.vdjca"), emit: vdjca, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def java_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ''
    """
    mixcr $java_mem analyze \\
        $preset \\
        --species $species \\
        $args \\
        --threads $task.cpus \\
        ${reads} \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr -v | sed -n '1p' | sed -E 's/MiXCR v([0-9\\.]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clones_TRA.tsv
    touch ${prefix}.clones_TRB.tsv
    touch ${prefix}.clns
    touch ${prefix}.vdjca
    touch ${prefix}.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixcr: \$(mixcr -v | sed -n '1p' | sed -E 's/MiXCR v([0-9\\.]+).*/\\1/')
    END_VERSIONS
    """
}
