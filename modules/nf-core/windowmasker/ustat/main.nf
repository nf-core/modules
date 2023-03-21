process WINDOWMASKER_USTAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0':
        'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    tuple val(meta), path(counts)
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("${output}")  , emit: intervals
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args         ?: ""
    def prefix  = task.ext.prefix       ?: "${meta.id}"
    def outfmt  = task.ext.suffix       ?: "interval"
    output  = "${prefix}.${outfmt}"

    """
    windowmasker -ustat \\
        ${counts} \\
        $args \\
        -in ${ref} \\
        -out ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1)
    END_VERSIONS
    """
}
