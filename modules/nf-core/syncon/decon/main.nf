process SYNCON_DECON {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/be85cb074d7868dc6fbfcdf0e6eabde839e61021665c8c2f856aa96dcbc3e568/data'
        : 'community.wave.seqera.io/library/syncon:1.0--454827d120667dec'}"

    input:
    tuple val(meta), path(asm_1khash), path(asm_1syncseq), path(reads_1khash)
    val write_profile

    output:
    tuple val(meta), path("*.decon.txt"), emit: syncon_decon_txt
    tuple val(meta), path("*.profile.txt"), emit: syncon_profile_txt
    tuple val("${task.process}"), val('syncon'), val("1.0"), topic: versions, emit: versions_syncon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def profile_arg = write_profile ? "-p ${prefix}.profile.txt" : ''
    """
    syncon decon \\
        ${args} \\
        -o ${prefix}.decon.txt \\
        ${profile_arg} \\
        ${asm_1khash} \\
        ${asm_1syncseq} \\
        ${reads_1khash}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def profile_output = write_profile ? "touch ${prefix}.profile.txt" : ''
    """
    echo ${args}

    touch ${prefix}.decon.txt
    ${profile_output}
    """
}
