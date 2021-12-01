process BCFTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bcftools=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0' :
        'quay.io/biocontainers/bcftools:1.13--h3a49de5_0' }"

    input:
    tuple val(meta), path(vcf)
    path fai
    path header

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def update_sequences = fai ? "-f $fai" : ""
    def new_header       = header ? "-h $header" : ""
    """
    bcftools \\
        reheader \\
        $update_sequences \\
        $new_header \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix}.vcf.gz \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
