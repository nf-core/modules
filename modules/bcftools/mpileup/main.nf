process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0' :
        'quay.io/biocontainers/bcftools:1.13--h3a49de5_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*.gz")      , emit: vcf
    tuple val(meta), path("*.tbi")     , emit: tbi
    tuple val(meta), path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "${meta.id}" > sample_name.list

    bcftools mpileup \\
        --fasta-ref $fasta \\
        $args \\
        $bam \\
        | bcftools call --output-type v $args2 \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z $args3

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
