process BISMARK_METHYLATIONEXTRACTOR {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data'
        : 'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bedGraph.gz"), emit: bedgraph
    tuple val(meta), path("*.txt.gz"), emit: methylation_calls
    tuple val(meta), path("*.cov.gz"), emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt"), emit: mbias
    tuple val("${task.process}"), val('bismark'), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), emit: versions_bismark, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    if (!args.contains('--multicore') && task.cpus >= 6) {
        args += " --multicore ${(task.cpus / 3) as int}"
    }
    // Only set buffer_size when there are more than 6.GB of memory available
    if (!args.contains('--buffer_size') && task.memory?.giga > 6) {
        args += " --buffer_size ${task.memory.giga - 2}G"
    }

    def seqtype = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        ${bam} \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --report \\
        ${seqtype} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bedGraph.gz
    echo "" | gzip > ${prefix}.txt.gz
    echo "" | gzip > ${prefix}.cov.gz
    touch ${prefix}_splitting_report.txt
    touch ${prefix}.M-bias.txt
    """
}
