process MTNUCRATIO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mtnucratio=0.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mtnucratio:0.7--hdfd78af_2"
    } else {
        container "quay.io/biocontainers/mtnucratio:0.7--hdfd78af_2"
    }

    input:
    tuple val(meta), path(bam)
    val(mt_id)

    output:
    tuple val(meta), path("*.mtnucratio"), emit: mtnucratio
    tuple val(meta), path("*.json")      , emit: json
    path "versions.yml"                  , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mtnucratio \\
        $args \\
        $bam \\
        $mt_id

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(mtnucratio --version 2>&1) | head -n1 | sed 's/Version: //')
    END_VERSIONS
    """
}
