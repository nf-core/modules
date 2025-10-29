process DEACON_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.5.0--h4349ce8_0':
        'biocontainers/deacon:0.5.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(index), path(fastq)

    output:
    tuple val(meta), path("*_filtered.fq"), emit: fastq_filtered
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    deacon \\
        filter \\
        --threads ${task.cpus} \\
        $args \\
        -d $index \\
        $fastq > ${prefix}_filtered.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fastq.baseName}"
    """
    touch ${prefix}_filtered.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """
}
