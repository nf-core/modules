process PASTY {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::pasty=1.0.0" : null)
    def container_image = "pasty:1.0.0--hdfd78af_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: tsv
    tuple val(meta), path("${prefix}.blastn.tsv") , emit: blast
    tuple val(meta), path("${prefix}.details.tsv"), emit: details
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pasty \\
        $args \\
        --prefix $prefix \\
        --assembly $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/^.*pasty, version //;' )
    END_VERSIONS
    """
}
