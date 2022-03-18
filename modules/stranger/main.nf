process STRANGER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::stranger=0.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stranger:0.8.1--pyh5e36f6f_0':
        'quay.io/biocontainers/stranger:0.8.1--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stranger \\
        $args \\
        $vcf | gzip --no-name > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stranger: \$( stranger --version )
    END_VERSIONS
    """
}
