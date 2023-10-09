process CNVNATOR_CONVERT2VCF {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::cnvnator=0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvnator:0.4.1--py310h2dce045_7':
        'biocontainers/cnvnator:0.4.1--py310h2dce045_7' }"

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
    def output_cmd = "> ${prefix}_cnvnator.vcf"
    """
    cnvnator2VCF.pl \\
        ${calls} \\
        $args \\
        $output_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
