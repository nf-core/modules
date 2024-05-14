process BLAST_TBLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0':
        'biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB=`find -L ./ -name "*.nsq" | sed 's/\\.nsq\$//'`
    tblastn \\
        -num_threads $task.cpus \\
        -db \$DB \\
        -query $fasta \\
        $args \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
