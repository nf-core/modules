process MODKIT_LOCALIZE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta), path(bedmethyl), path(tbi)
    tuple val(meta2), path(sizes)
    tuple val(meta3), path(regions)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv  , optional: true
    tuple val(meta), path("*.html") , emit: chart , optional: true
    tuple val(meta), path("*.log")  , emit: log   , optional: true
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def chart   = "--chart ${prefix}.localize.html --name ${prefix}"

    """
    modkit \\
        localize \\
        $args \\
        --threads ${task.cpus} \\
        --genome-sizes $sizes \\
        --regions $regions \\
        $chart \\
        --force \\
        -o ${prefix}.localize.tsv \\
        $bedmethyl
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.localize.tsv
    touch ${prefix}.localize.html
    touch ${prefix}.log
    """
}
