process SEQKIT_SLIDING {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.fast*"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta$|.+\.fa$|.+\.fas$|.+\.fna$/) {
        extension = "fasta"
    }
    """
    seqkit \\
        sliding \\
        ${fastx} \\
        ${args} \\
        --threads ${task.cpus} \\
        -o ${prefix}.${extension}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fastx}" ==~ /.+\.fasta$|.+\.fa$|.+\.fas$|.+\.fna$/) {
        extension = "fasta"
    }
    """
    touch ${prefix}.${extension}
    """
}
