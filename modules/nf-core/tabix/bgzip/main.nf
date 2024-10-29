process TABIX_BGZIP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2' :
        'biocontainers/htslib:1.20--h5efdd21_2' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${output}")    , emit: output
    tuple val(meta), path("${output}.gzi"), emit: gzi, optional: true
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    in_bgzip = ["gz", "bgz", "bgzf"].contains(input.getExtension())
    extension = in_bgzip ? input.getBaseName().tokenize(".")[-1] : input.getExtension()
    output   = in_bgzip ? "${prefix}.${extension}" : "${prefix}.${extension}.gz"
    command = in_bgzip ? '-d' : ''
    // Name the index according to $prefix, unless a name has been requested
    if ((args.matches("(^| )-i\\b") || args.matches("(^| )--index(\$| )")) && !args.matches("(^| )-I\\b") && !args.matches("(^| )--index-name\\b")) {
        args = args + " -I ${output}.gzi"
    }
    """
    bgzip $command -c $args -@${task.cpus} $input > ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    in_bgzip = ["gz", "bgz", "bgzf"].contains(input.getExtension())
    output   = in_bgzip ? input.getBaseName() : "${prefix}.${input.getExtension()}.gz"

    """
    echo "" | gzip > ${output}
    touch ${output}.gzi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
