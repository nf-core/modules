process ASSEMBLYSCAN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::assembly-scan=0.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/assembly-scan:0.4.1--pyhdfd78af_0' :
        'quay.io/biocontainers/assembly-scan:0.4.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assembly-scan $assembly > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemblyscan: \$( assembly-scan --version 2>&1 | sed 's/^.*assembly-scan //; s/Using.*\$//' )
    END_VERSIONS
    """
}
