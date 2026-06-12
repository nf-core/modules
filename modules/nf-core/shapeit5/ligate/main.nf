process SHAPEIT5_LIGATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa06870e893f9045944461e5b674adffec6deec41e8496b63ec54d44a98d6134/data'
        : 'community.wave.seqera.io/library/shapeit5:5.1.1--09a6cb254ece8f6e'}"

    input:
    tuple val(meta), path(input_list), path(input_list_index)
    val suffix_

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    tuple val("${task.process}"), val('shapeit5'), eval('SHAPEIT5_ligate | sed "5!d;s/^.*Version *: //; s/ .*$//"'), topic: versions, emit: versions_shapeit5


    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = suffix_ ?: "vcf.gz"
    assert suffix in ["vcf", "vcf.gz", "bcf", "bcf.gz"]

    // Create file list using Groovy (most portable)
    def file_list = input_list
        .collect { file_path -> file_path.toString() }
        .sort()
        .join('\n')

    """
    echo "${file_list}" > all_files.txt

    SHAPEIT5_ligate \\
        ${args} \\
        --input all_files.txt \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = suffix_ ?: "vcf.gz"
    assert suffix in ["vcf", "vcf.gz", "bcf", "bcf.gz"]
    def create_cmd = suffix.endsWith(".gz") ? 'echo "" | gzip >' : "touch"
    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
