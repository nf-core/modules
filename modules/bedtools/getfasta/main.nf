process BEDTOOLS_GETFASTA {
    tag "$bed"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    path bed
    path fasta

    output:
    path "*.fa"         , emit: fasta
    path "versions.yml" , emit: versions

    script:
    def prefix   = options.suffix ? "${bed.baseName}${options.suffix}" : "${bed.baseName}"
    """
    bedtools \\
        getfasta \\
        $options.args \\
        -fi $fasta \\
        -bed $bed \\
        -fo ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
