process GRAPHTYPER_VCFCONCATENATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::graphtyper=2.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphtyper:2.7.2--h7d7f7ad_0':
        'quay.io/biocontainers/graphtyper:2.7.2--h7d7f7ad_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.concatenated.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graphtyper vcf_concatenate \\
        $vcf \\
        $args \\
        --write_tbi \\
        --output=${prefix}.concatenated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
    END_VERSIONS
    """
}
