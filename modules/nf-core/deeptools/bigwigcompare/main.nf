process DEEPTOOLS_BIGWIGCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(bigwig1)     , path(bigwig2)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.{bigWig,bedgraph}"), emit: output
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args                                  ?: ""
    def prefix = task.ext.prefix                              ?: "${meta.id}"
    def blacklist_cmd = blacklist                             ? "--blackListFileName ${blacklist}" : ""
    def extension = args.contains("--outFileFormat bedgraph") ? "bedgraph"                         : "bigWig"

    """
    bigwigCompare \\
        --bigwig1 $bigwig1 \\
        --bigwig2 $bigwig2 \\
        --outFileName ${prefix}.${extension} \\
        --numberOfProcessors $task.cpus \\
        $blacklist_cmd \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigCompare --version | sed -e "s/bigwigCompare //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args                                  ?: ''
    def prefix = task.ext.prefix                              ?: "${meta.id}"
    def extension = args.contains("--outFileFormat bedgraph") ? "bedgraph" : "bigWig"

    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigCompare --version | sed -e "s/bigwigCompare //g")
    END_VERSIONS
    """
}
