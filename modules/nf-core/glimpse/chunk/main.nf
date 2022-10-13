process GLIMPSE_CHUNK {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::glimpse-bio=1.1.1" : null)
    def container_image = "/glimpse-bio:1.1.1--h2ce4488_2"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(input)
    val region

    output:
    tuple val(meta), path("*.txt"), emit: chunk_chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ""
    def VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    GLIMPSE_chunk \\
        $args \\
        --input $input \\
        --region $region \\
        --thread $task.cpus \\
        --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: $VERSION
    END_VERSIONS
    """
}
