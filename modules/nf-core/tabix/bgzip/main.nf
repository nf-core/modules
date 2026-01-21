process TABIX_BGZIP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${output}"), emit: output
    tuple val(meta), path("*.gzi")    , emit: gzi, optional: true
    path  "versions.yml"              , emit: versions

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
    split_args = args.split(' +|=')
    if ((split_args.contains('-i') || split_args.contains('--index')) && !split_args.contains('-I') && !split_args.contains('--index-name')) {
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
