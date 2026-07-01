process CATPACK_SUMMARISE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e58abd32df3c4d65314311ad1b3b8355fa3d9f229aa3e39a1db6c938eb406d1/data'
        : 'community.wave.seqera.io/library/cat:6.0.1--6bec964806e8076a'}"

    input:
    tuple val(meta), path(classification)
    tuple val(meta2), path(contigs)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('catpack'), eval("CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g'"), topic: versions, emit: versions_catpack

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${classification}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    def insert_contigs = contigs ? "-c ${contigs}" : ''
    """
    CAT_pack summarise \\
        ${args} \\
        -i ${classification} \\
        ${insert_contigs} \\
        -o ${prefix}.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${classification}" == "${prefix}.txt") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    def insert_contigs = contigs ? "-c ${contigs}" : ''
    """
    echo "CAT_pack summarise \\
        ${args} \\
        -i ${classification} \\
        ${insert_contigs} \\
        -o ${prefix}.txt"

    touch ${prefix}.txt
    """
}
