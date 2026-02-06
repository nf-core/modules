process GLIMPSE_SAMPLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--hce55b13_1'
        : 'biocontainers/glimpse-bio:1.1.1--hce55b13_1'}"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: haplo_sampled
    tuple val("${task.process}"), val('glimpse'), eval("GLIMPSE_sample --help | sed -n '/Version/s/.*: //p'"), topic: versions, emit: versions_glimpse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    """
    GLIMPSE_sample \\
        ${args} \\
        --input ${input} \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    def create_cmd = suffix.endsWith(".gz") ? "echo | gzip >" : "touch"

    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
