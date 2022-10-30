process COOLER_BALANCE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0':
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(cool), val(resolution)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: cool
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = resolution ? "::/resolutions/$resolution" : ""
    extension = cool.getExtension()
    if ("$cool" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    cp ${cool} ${prefix}.${extension}

    cooler balance \\
        $args \\
        -p ${task.cpus} \\
        ${prefix}.${extension}${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
