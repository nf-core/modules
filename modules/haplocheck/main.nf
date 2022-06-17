process HAPLOCHECK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::haplocheck=1.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplocheck:1.3.3--h4a94de4_0':
        'quay.io/biocontainers/haplocheck:1.3.3--h4a94de4_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.txt") , emit: txt
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplocheck --raw --out $prefix $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplocheck: \$(echo \$(haplocheck --version 2>&1) | cut -f 2 -d " " )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.txt
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplocheck: \$(echo \$(haplocheck --version 2>&1) | cut -f 2 -d " " )
    END_VERSIONS
    """
}
