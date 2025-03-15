process STIMULUS_TRANSFORM {
    tag "$meta.id"
    label 'process_medium'
    container "docker.io/luisas/stimulus-py:0.0.10"

    input:
    tuple val(meta), path(transform_json)
    path(csv)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: csv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def STIMULUS_VER = '0.0.9'
    """
    stimulus-transform-csv -c ${csv} -j ${transform_json} -o ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$( pip show stimulus-py | grep Version | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def STIMULUS_VER = '0.0.9' // container not used in stub, change manually
    """
    touch ${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
