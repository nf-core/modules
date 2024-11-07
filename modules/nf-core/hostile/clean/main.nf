process HOSTILE_CLEAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:1.1.0--pyhdfd78af_0':
        'biocontainers/hostile:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("cleaned_reads/"), emit: cleaned_reads
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_cmd = meta.single_end ? "--fastq1 ${reads.sort()[0]}" : "--fastq1 ${reads.sort()[0]} --fastq2 ${reads.sort()[1]}"
    """
    export HOSTILE_CACHE_DIR=${reference}
    mkdir cleaned_reads/

    hostile \\
        clean \\
        $args \\
        --threads $task.cpus \\
        --out-dir cleaned_reads/ \\
        --reorder \\ # For reproducability
        --offline \\ # This module requires input reference files so never download in this process
        $reads_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir cleaned_reads/
    touch ${prefix}_r1.fastq.gz
    touch ${prefix}_r2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """
}
