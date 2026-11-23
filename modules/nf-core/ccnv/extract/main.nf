process CCNV_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd552e195217d37801cf86348dbe875fa66ae2014147e481fd6be6e80dbbcd90/data'
        : 'community.wave.seqera.io/library/ccnv:1.1.0--726d8b142feb0e3d'}"

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

    if (!reads_index) {
        error "CCNV_EXTRACT: Missing reads index file"
    }
    """
    ccnv extract\\
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
