process FGUMI_EXTRACT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a510706f3481fae12ff6100d6e4ad298b8bf464a2d93a6afe35e9cf26542d080/data'
        : 'community.wave.seqera.io/library/fgumi:0.2.0--fe028e7a64e5da27'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!args.contains('--sample') || !args.contains('--library')) {
        error("fgumi extract requires both --sample and --library to be supplied via task.ext.args")
    }

    """
    fgumi extract \\
        --inputs ${reads.join(' ')} \\
        --output ${prefix}.bam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
