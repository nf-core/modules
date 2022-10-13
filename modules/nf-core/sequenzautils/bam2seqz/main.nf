process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    def container_image = "/sequenza-utils:3.0.0--py38h6ed170a_2"
                                                         container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(normalbam), path(tumourbam)
    path fasta
    path wigfile

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        bam2seqz \\
        $args \\
        -n $normalbam \\
        -t $tumourbam \\
        --fasta $fasta \\
        -gc $wigfile \\
        -o ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
