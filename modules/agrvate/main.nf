process AGRVATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::agrvate=1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/agrvate:1.0.2--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/agrvate:1.0.2--hdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}-results/${fasta.baseName}-summary.tab"), emit: summary
    path "${fasta.baseName}-results"                                                , emit: results_dir
    path "versions.yml"                                                             , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    agrvate \\
        $options.args \\
        -i $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
