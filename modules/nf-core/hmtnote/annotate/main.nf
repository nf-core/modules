process HMTNOTE_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmtnote:0.7.2--pyhdfd78af_1':
        'biocontainers/hmtnote:0.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    tuple val("${task.process}"), val("hmtnote"), eval("hmtnote --version | sed 's/^.*hmtnote, version //'"), emit: versions_hmtnote, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    hmtnote \\
        annotate \\
        $vcf \\
        ${prefix}.vcf \\
        $args
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
