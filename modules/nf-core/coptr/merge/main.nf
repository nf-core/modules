process COPTR_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coptr:1.1.4--pyhdfd78af_3':
        'biocontainers/coptr:1.1.4--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def input_bams  = "${bams.join(" ")}"
    """
    coptr \
        merge \
        $args \
        ${input_bams} \
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """
}
