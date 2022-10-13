process ARTIC_GUPPYPLEX {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::artic=1.2.2" : null)
    def container_image = "/artic:1.2.2--pyhdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.2' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
    """
    artic \\
        guppyplex \\
        $args \\
        --directory $fastq_dir \\
        --output ${prefix}.fastq

    pigz -p $task.cpus *.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: $VERSION
    END_VERSIONS
    """
}
