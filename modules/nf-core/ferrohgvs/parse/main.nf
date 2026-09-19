process FERROHGVS_PARSE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e0f241c7350e07e9ed616a5402ce084099a3614dd26b8afee6a02aa041e5fd5/data'
        : 'community.wave.seqera.io/library/ferro-hgvs:0.6.0--80a6f47dc50a3cf3'}"

    input:
    tuple val(meta), path(variants)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: parsed
    tuple val("${task.process}"), val('ferrohgvs'), eval("ferro --version | sed 's/^ferro //'"), emit: versions_ferrohgvs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = args.contains("--format json") || args.contains("-f json") ? "json" : "txt"
    if ("${variants}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    ferro \\
        parse \\
        --input ${variants} \\
        --output ${prefix}.${suffix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = args.contains("--format json") || args.contains("-f json") ? "json" : "txt"
    """
    touch ${prefix}.${suffix}
    """
}
