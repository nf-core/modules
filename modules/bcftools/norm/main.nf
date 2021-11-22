process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
    }

    input:
    tuple val(meta), path(vcf)
    path(fasta)

    output:
    tuple val(meta), path("*.gz") , emit: vcf
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --output ${prefix}.vcf.gz \\
        $options.args \\
        --threads $task.cpus \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
