process CCNV_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f3fc624233aa2648bc513082af307dd7061f5fd25f4cd53fc749f8c62050496b/data'
        : 'community.wave.seqera.io/library/ccnv:1.0.0--c01430de16e6765e'}"

    input:
    tuple val(meta), path(reads), path(reads_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.ccnv"), emit: ccnv
    tuple val("${task.process}"), val('ccnv'), eval("ccnv --version 2>&1 | sed 's/ccnv //'"), topic: versions, emit: versions_ccnv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta ? "--fasta $fasta" : ""
    """
    ccnv \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        $fasta_arg \\
        $args \\
        $reads
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.ccnv
    """
}
