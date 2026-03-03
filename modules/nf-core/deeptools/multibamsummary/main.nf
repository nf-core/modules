process DEEPTOOLS_MULTIBAMSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(bams)    , path(bais), val(labels)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.npz"), emit: matrix
    tuple val("${task.process}"), val('deeptools'), eval('multiBamSummary --version | sed "s/multiBamSummary //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "all_bam"
    def blacklist_cmd = blacklist ? "--blackListFileName ${blacklist}" : ""
    def label  = labels ? "--labels ${labels.join(' ')}" : ''
    """
    multiBamSummary bins \\
        $args \\
        $label \\
        --bamfiles ${bams.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}.bamSummary.npz \\
        $blacklist_cmd
    """

    stub:
    def prefix = task.ext.prefix ?: "all_bam"
    """
    touch ${prefix}.bamSummary.npz
    """
}
