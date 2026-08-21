process SEMIBIN_CONCATENATEFASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4c049230e87478349b2ecca53156faedfb077df85a4f6deaa8d0b22bd8a975c4/data'
        : 'community.wave.seqera.io/library/semibin:2.4.1--594344a862aee1b8'}"

    input:
    tuple val(meta), path(fastas)
    val compress_out

    output:
    tuple val(meta), path("*.fasta*"), emit: concat_fasta
    tuple val("${task.process}"), val('semibin'), eval("SemiBin2 --version"), topic: versions, emit: versions_semibin

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile = compress_out ? "${prefix}.fasta.gz" : "${prefix}.fasta"
    """
    SemiBin2 \\
        ${args} \\
        concatenate_fasta \\
        ${args2} \\
        --input-fasta ${fastas} \\
        --output ${outfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stub_command = compress_out ? "echo '' | gzip > ${prefix}.fasta.gz" : "echo '' > ${prefix}.fasta"
    """
    ${stub_command}
    """
}
