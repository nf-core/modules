
process SYLPHTAX_TAXPROF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.2.0--pyhdfd78af_0':
        'biocontainers/sylph-tax:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sylph_results)
    path taxonomy

    output:
    tuple val(meta), path("*.sylphmpa"), emit: taxprof_output
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    sylph-tax \\
        taxprof \\
        $sylph_results \\
        $args \\
        -t $taxonomy

    mv *.sylphmpa ${prefix}.sylphmpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    touch ${prefix}.sylphmpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version)
    END_VERSIONS
    """
}
