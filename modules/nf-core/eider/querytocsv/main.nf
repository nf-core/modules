process EIDER_QUERYTOCSV {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eider:0.3--hdfd78af_0' :
        'quay.io/biocontainers/eider:0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(sql)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('eider'), eval("eider --version 2>&1 | grep -o 'eider .*' | cut -f2 -d ' '"), topic: versions, emit: versions_eider

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    eider \
        $args \
        --skip-history \
        --verbose \
        --parameters prefix=$prefix \
        --query-path $sql
    """
}
