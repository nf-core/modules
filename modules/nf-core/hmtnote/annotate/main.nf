process HMTNOTE_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::hmtnote=0.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmtnote:0.7.2--pyhdfd78af_1':
        'biocontainers/hmtnote:0.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_annotated.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    hmtnote \\
        annotate \\
        $vcf \\
        ${prefix}_annotated.vcf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmtnote: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' ))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmtnote: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
