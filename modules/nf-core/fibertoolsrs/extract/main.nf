process FIBERTOOLSRS_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fibertools-rs:0.6.2--h3b373d1_0':
        'biocontainers/fibertools-rs:0.7.1--h3b373d1_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outbed = params.ft_extract_type ? "--${params.ft_extract_type} ${prefix}.bed" : "--all ${prefix}.bed"

    """
    ft \\
        extract \\
        $args \\
        --threads $task.cpus \\
        $bam \\
        $outbed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fibertools-rs: \$(ft --version | sed -E 's/.* ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fibertools-rs: \$(ft --version | sed -E 's/.* ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/' )
    END_VERSIONS
    """
}
