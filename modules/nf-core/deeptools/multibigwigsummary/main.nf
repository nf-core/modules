process DEEPTOOLS_MULTIBIGWIGSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba9641d61acf134c6b989b3af5f25bcd5aadb3faf692a41c00ba4cd25f99032e/data'
:         'community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3' }"

    input:
    tuple val(meta) , path(bigwigs) , val(labels)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.npz"), emit: matrix
    tuple val("${task.process}"), val('deeptools'), eval('multiBigwigSummary --version | sed "s/multiBigwigSummary //g"') , emit: versions_deeptools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "all_bigwig"
    """
    touch ${prefix}.bigwigSummary.npz
    """
}
