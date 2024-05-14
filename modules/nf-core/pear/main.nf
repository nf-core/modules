process PEAR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pear=0.9.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pear:0.9.6--h67092d7_8':
        'biocontainers/pear:0.9.6--h67092d7_8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.assembled.fastq.gz")                                                  , emit: assembled
    tuple val(meta), path("*.unassembled.forward.fastq.gz"), path("*.unassembled.reverse.fastq.gz"), emit: unassembled
    tuple val(meta), path("*.discarded.fastq.gz")                                                  , emit: discarded
    path "versions.yml"                                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -f ${reads[0]}
    gunzip -f ${reads[1]}
    pear \\
        -f ${reads[0].baseName} \\
        -r ${reads[1].baseName} \\
        -o $prefix \\
        -j $task.cpus \\
        $args
    gzip -f ${prefix}.assembled.fastq
    gzip -f ${prefix}.unassembled.forward.fastq
    gzip -f ${prefix}.unassembled.reverse.fastq
    gzip -f ${prefix}.discarded.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pear: \$(pear -h | grep 'PEAR v' | sed 's/PEAR v//' | sed 's/ .*//' ))
    END_VERSIONS
    """
}
