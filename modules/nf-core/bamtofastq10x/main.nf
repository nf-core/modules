process BAMTOFASTQ10X {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/10x_bamtofastq:1.4.1--hdbdd923_2':
        'biocontainers/10x_bamtofastq:1.4.1--hdbdd923_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("**/*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('bamtofastq10x'), eval('bamtofastq --version |& sed "1!d ; s/bamtofastq //"'), emit: versions_bamtofastq10x, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamtofastq \\
        $args \\
        $bam \\
        $prefix
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/bamtofastq10x
    echo "" | gzip > ${prefix}/bamtofastq10x/bamtofastq.fastq.gz
    """
}
