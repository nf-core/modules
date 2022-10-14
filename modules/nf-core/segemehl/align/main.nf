process SEGEMEHL_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::segemehl=0.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5':
        'quay.io/biocontainers/segemehl:0.3.4--hc2ea5fd_5' }"

    input:
    tuple val(meta), path(reads)
    path(fasta)
    path(index)

    output:
    tuple val(meta), path("${meta.id}/*"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads = meta.single_end ? "-q ${reads}" : "-q ${reads[0]} -p ${reads[1]}"
    """
    mkdir -p $prefix

    segemehl.x \\
        -t $task.cpus \\
        -d $fasta \\
        -i $index \\
        $reads \\
        $args \\
        -o ${prefix}/${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    touch $prefix/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """
}
