// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GLIMPSE_LIGATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f722eb1dd56d995b28e1ec47f7808450bf6477e0b4ec28ec5f294c1a66fd479e/data'
        : 'community.wave.seqera.io/library/glimpse-bio:1.1.1--26d63752054f98c5'}"

    input:
    tuple val(meta), path(input_list), path(input_index)
    val extension

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    tuple val("${task.process}"), val('glimpse'), eval("GLIMPSE_ligate --help | sed -n '/Version/s/.*: //p'"), topic: versions, emit: versions_glimpse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = extension ?: "vcf.gz"

    // Create file list using Groovy (most portable)
    def file_list = input_list
        .collect { file_path -> file_path.toString() }
        .sort()
        .join('\n')

    """
    echo "${file_list}" > all_files.txt

    GLIMPSE_ligate \\
        ${args} \\
        --input all_files.txt \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = extension ?: "vcf.gz"

    def create_cmd = suffix.endsWith(".gz") ? "echo '' | gzip >" : "touch"

    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
