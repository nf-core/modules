process STIMULUS_SPLITCSV {
    tag "$meta.id"
    label 'process_low'

    // #TO-DO: UPDATE
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    conda "${moduleDir}/environment.yml"
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    tuple val(meta), path(split_json)
    path(data_csv)

    output:
    tuple val(meta), path("${meta.id}.csv"), emit: split_csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    stimulus-split-csv -c ${data_csv} -j ${split_json} -o ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$( pip show stimulus-py | grep Version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def STIMULUS_VER = '0.0.9' // container not used in stub, change manually
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
