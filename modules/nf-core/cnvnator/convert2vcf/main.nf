process CNVNATOR_CONVERT2VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvnator:0.4.1--py312h99c8fb2_11':
        'biocontainers/cnvnator:0.4.1--py312h99c8fb2_11' }"

    input:
    tuple val(meta), path(calls)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    cnvnator2VCF.pl \\
        ${calls} \\
        $args \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNVnator : \$(echo \$(cnvnator 2>&1 | sed -n '3p' | sed 's/CNVnator v//'))
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNVnator : \$(echo \$(cnvnator 2>&1 | sed -n '3p' | sed 's/CNVnator v//'))
    END_VERSIONS
    """
}
