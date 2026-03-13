process SKESA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa%3A2.5.1--h077b44d_3':
        'biocontainers/skesa:2.5.1--h077b44d_3' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('skesa'), eval("skesa --version 2>&1 | sed -n 's/^SKESA //p'"), topic: versions, emit: versions_skesa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory ? "--memory ${task.memory.toMega()} MB" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? "" : "--use_paired_ends"
    """
        skesa \
            --reads $fastq \
            --contigs_out ${prefix}.fa \
            --cores $task.cpus \
            $paired_end \
            $memory \
            $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}.fa
    """
}
