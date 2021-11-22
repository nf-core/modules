process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    samtools \\
        depth \\
        $args \\
        -o ${prefix}.tsv \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
