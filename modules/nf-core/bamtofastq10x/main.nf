process BAMTOFASTQ10X {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95e369e8e54b4d483c85a1584cf0874f5d2efd547928223557e3d8da2f7ee5ec/data'
        : 'community.wave.seqera.io/library/10x_bamtofastq:1.4.1--a7cbf7b60c2c3fd7' }"

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
