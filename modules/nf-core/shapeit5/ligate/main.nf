process SHAPEIT5_LIGATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0'
        : 'biocontainers/shapeit5:5.1.1--hb60d31d_0'}"

    input:
    tuple val(meta), path(input_list), path(input_list_index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    tuple val("${task.process}"), val('shapeit5'), eval('SHAPEIT5_ligate | sed "5!d;s/^.*Version *: //; s/ .*$//"'), topic: versions, emit: versions_shapeit5

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

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
    def suffix = task.ext.suffix ?: "vcf.gz"

    def create_cmd = suffix.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
