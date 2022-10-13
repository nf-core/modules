process SEQUENZAUTILS_GCWIGGLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
        def container_image = "/sequenza-utils:3.0.0--py38h6ed170a_2"
                                                         container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig.gz"), emit: wig
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        gc_wiggle \\
        $args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
