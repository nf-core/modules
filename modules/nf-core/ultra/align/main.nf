process ULTRA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ultra_bioinformatics=0.1 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0':
        'biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6:c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0' }"

    input:
    tuple val(meta), path(reads)
    path genome
    tuple path(pickle), path(db)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
