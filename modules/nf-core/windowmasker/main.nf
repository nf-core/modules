process WINDOWMASKER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0':
        'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*")    , emit: wm_intervals
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1   = task.ext.args1    ?: ""
    def args2   = task.ext.args2    ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def outfmt  = task.ext.outfmt   ?: "interval"

    """
    windowmasker -mk_counts \\
        $args1 \\
        -in ${ref} \\
        -out ${prefix}.txt

    windowmasker -ustat \\
        ${prefix}.txt \\
        $args2 \\
        -in ${ref} \\
        -out ${prefix}.${outfmt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full 2>&1)
    END_VERSIONS
    """
}
