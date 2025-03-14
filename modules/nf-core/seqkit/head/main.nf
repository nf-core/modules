process SEQKIT_HEAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.0--h9ee0642_0':
        'biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastqs), val(seq_count)

    output:
    tuple val(meta), path("*_subset.fastq.gz"), emit: subset
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    for f in ${fastqs.join(' ')}
    do
        seqkit head \\
            ${args} \\
            --threads $task.cpus \\
            -n ${seq_count} \\
            -o "\$(basename \$f .fastq.gz)_subset.fastq.gz" \\
            \$f
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    """
    for f in ${fastqs.join(' ')}
    do
       echo "" | gzip > "\$(basename \$f .fastq.gz)_subset.fastq.gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
