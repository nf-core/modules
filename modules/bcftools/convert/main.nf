process BCFTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0':
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(bed)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def regions = bed ? "--regions-file $bed" : ""
    def reference = fasta ?  "--fasta-ref $fasta" : ""

    """
    bcftools convert \\
        $args \\
        $regions \\
        --output ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        $reference \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
