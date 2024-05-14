process BCFTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: split_vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools index -s ${vcf} |\\
    cut -f 1 | \\
    while read C; do \\
        bcftools view -O z -o ${prefix}.\${C}.vcf.gz ${vcf} --regions "\${C}" ;\\
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch "${prefix}.chr1.vcf.gz"
    touch "${prefix}.chr2.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
