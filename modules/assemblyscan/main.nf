process ASSEMBLYSCAN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::assembly-scan=0.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/assembly-scan:0.4.1--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/assembly-scan:0.4.1--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    assembly-scan $assembly > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( assembly-scan --version 2>&1 | sed 's/^.*assembly-scan //; s/Using.*\$//' )
    END_VERSIONS
    """
}
