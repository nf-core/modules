process LDAK_FILTERRELATEDNESS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(grm_details), path(grm_adjust)

    output:
    tuple val(meta), path("${prefix}.keep"), path("${prefix}.lose"), path("${prefix}.maxrel"), emit: filtered_list
    tuple val(meta), path("${prefix}_unrel.grm.bin"), path("${prefix}_unrel.grm.id"), path("${prefix}_unrel.grm.details"), path("${prefix}_unrel.grm.adjust"), emit: filtered_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def grm_prefix = grm_bin.name.replaceFirst(/\.grm\.bin$/, '')
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    ldak6 \
        --filter ${prefix} \
        --grm ${grm_prefix} \
        --max-threads ${task.cpus} \
        ${args}

    ldak6 \
        --sub-grm ${prefix}_unrel \
        --grm ${grm_prefix} \
        --keep ${prefix}.keep \
        --max-threads ${task.cpus} \
        ${args2}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.keep
    touch ${prefix}.lose
    echo "0.0" > ${prefix}.maxrel
    touch ${prefix}_unrel.grm.bin
    touch ${prefix}_unrel.grm.id
    touch ${prefix}_unrel.grm.details
    touch ${prefix}_unrel.grm.adjust
    """
}
