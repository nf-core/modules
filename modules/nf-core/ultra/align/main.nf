process ULTRA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0':
        'quay.io/biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(genome)
    tuple val(meta3), path(pickle), path(db)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('ultra'), eval("uLTRA --version | sed 's/uLTRA //'"), emit: versions_ultra, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    uLTRA \\
        align \\
        --t $task.cpus \\
        --prefix $prefix \\
        --index ./ \\
        $args \\
        $genome \\
        $reads \\
        ./

    samtools \\
        sort \\
        --threads $task.cpus \\
        -o ${prefix}.bam \\
        -O BAM \\
        $args2 \\
        ${prefix}.sam

    rm ${prefix}.sam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
