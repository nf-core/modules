process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: metrics
    tuple val(meta), path("*.pdf"), emit: histogram
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectInsertSizeMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectInsertSizeMetrics \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.txt \\
        --Histogram_FILE ${prefix}.pdf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        log.info '[Picard CollectInsertSizeMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    """
    touch ${prefix}.pdf
    touch ${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """



}
