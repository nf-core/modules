process MENINGOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::meningotype=0.8.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/meningotype:0.8.5--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/meningotype:0.8.5--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    meningotype \\
        $options.args \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(meningotype --version 2>&1) | sed 's/^.*meningotype v//' )
    END_VERSIONS
    """
}
