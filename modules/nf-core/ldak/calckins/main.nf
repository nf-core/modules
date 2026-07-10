process LDAK_CALCKINS {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(weights_file)
    val power

    output:
    tuple val(meta), path("${prefix}.grm.bin"), path("${prefix}.grm.id"), path("${prefix}.grm.details"), path("${prefix}.grm.adjust"), emit: ldak_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def weights_arg = weights_file ? "--weights ${weights_file}" : "--ignore-weights YES"
    """
    ldak6 \\
        --calc-kins-direct ${prefix} \\
        --bfile ${bed.baseName} \\
        --power ${power} \\
        ${weights_arg} \\
        --max-threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.id
    touch ${prefix}.grm.details
    touch ${prefix}.grm.adjust
    """
}
