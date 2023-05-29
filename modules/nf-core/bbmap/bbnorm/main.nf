process BBMAP_BBNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01 pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0':
        'biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input  = meta.single_end ? "in=${fastq.join(',')}" : "in=${fastq[0]} in2=${fastq[1]}"
    output = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.nm.fastq.gz out2=${prefix}_2.nm.fastq.gz"

    """
    bbnorm.sh \\
        $input \\
        $output \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbnorm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
