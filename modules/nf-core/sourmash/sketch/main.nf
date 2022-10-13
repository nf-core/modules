process SOURMASH_SKETCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::sourmash=4.2.4" : null)
    def container_image = "/sourmash:4.2.4--hdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "dna --param-string 'scaled=1000,k=31'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sketch \\
        $args \\
        --merge '${prefix}' \\
        --output '${prefix}.sig' \\
        $sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
