process LISSERO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::lissero=0.4.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lissero:0.4.9--py_0' :
        'quay.io/biocontainers/lissero:0.4.9--py_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lissero \\
        $args \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lissero: \$( echo \$(lissero --version 2>&1) | sed 's/^.*LisSero //' )
    END_VERSIONS
    """
}
