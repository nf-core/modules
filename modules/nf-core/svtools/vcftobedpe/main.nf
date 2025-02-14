process SVTOOLS_VCFTOBEDPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtools:0.5.1--py_0':
        'biocontainers/svtools:0.5.1--py_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    svtools vcftobedpe \\
        $args
        --input $vcf \\
        --output ${prefix}.bedpe \\
        --tempdir ./tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtools: \$(svtools --version |& sed '1!d ; s/svtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedpe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtools: \$(svtools --version |& sed '1!d ; s/svtools //')
    END_VERSIONS
    """
}
