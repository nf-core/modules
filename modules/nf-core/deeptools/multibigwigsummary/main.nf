process DEEPTOOLS_MULTIBIGWIGSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(bigwigs) , val(labels)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.npz"), emit: matrix
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "all_bigwig"
    def blacklist_cmd = blacklist ? "--blackListFileName ${blacklist}" : ""
    def label  = labels ? "--labels ${labels.join(' ')}" : ''
    """
    multiBigwigSummary bins \\
        $args \\
        $label \\
        --bwfiles ${bigwigs.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}.bigwigSummary.npz \\
        $blacklist_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed -e "s/multiBigwigSummary //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_bigwig"
    """
    touch ${prefix}.bigwigSummary.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed -e "s/multiBigwigSummary //g")
    END_VERSIONS
    """
}
