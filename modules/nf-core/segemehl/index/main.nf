process SEGEMEHL_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::segemehl=0.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5':
        'biocontainers/segemehl:0.3.4--hc2ea5fd_5' }"

    input:
    path fasta

    output:
    path "*.idx",        emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${fasta.baseName}"
    """
    segemehl.x \\
        -t $task.cpus \\
        -d $fasta \\
        -x ${prefix}.idx \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """

    stub:
    def prefix = "${fasta.baseName}"
    """
    touch ${prefix}.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segemehl: \$(echo \$(segemehl.x 2>&1 | grep "ge5dee" | awk -F Z '{print substr(\$1, 2, 6)}' ))
    END_VERSIONS
    """
}
