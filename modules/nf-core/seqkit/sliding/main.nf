process SEQKIT_SLIDING {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.fast*"), emit: fastx
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta$|.+\.fa$|.+\.fas$|.+\.fna$/) {
        extension = "fasta"
    }
    """
    seqkit \\
        sliding \\
        ${fastx} \\
        ${args} \\
        --threads ${task.cpus} \\
        -o ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fastx" ==~ /.+\.fasta$|.+\.fa$|.+\.fas$|.+\.fna$/) {
        extension = "fasta"
    }
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
