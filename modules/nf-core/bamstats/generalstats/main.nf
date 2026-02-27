process BAMSTATS_GENERALSTATS {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bamstats:0.3.5--he881be0_0'
        : 'biocontainers/bamstats:0.3.5--he881be0_0'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val("bamstats"), eval('bamstats --version | grep "version: " | sed -e s"/version: //"'), topic: versions, emit: versions_bamstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Several ARGS are available.
    //  -a is a helpfu; one where you can add a BED file
    //  -u, --uniq outputs genomic coverage statistics for uniqely mapped reads
    """
    bamstats \\
        -i ${bam} \\
        ${args} \\
        -c ${task.cpus} \\
        -o ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.json
    """
}
