process SEGEMEHL_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::segemehl=0.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5':
        'biocontainers/segemehl:0.3.4--hc2ea5fd_5' }"

    input:
    tuple val(meta), path(reads)
    path(fasta)
    path(index)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads = meta.single_end ? "-q ${reads}" : "-q ${reads[0]} -p ${reads[1]}"
    def suffix = ( args.contains("-b") || args.contains("--bamabafixoida") ) ? "bam" : "sam"
    """
    mkdir -p $prefix

    segemehl.x \\
        -t $task.cpus \\
        -d $fasta \\
        -i $index \\
        $reads \\
        $args \\
        -o ${prefix}/${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = ( args.contains("-b") || args.contains("--bamabafixoida") ) ? "bam" : "sam"
    """
    mkdir -p $prefix
    touch ${prefix}/${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """
}
