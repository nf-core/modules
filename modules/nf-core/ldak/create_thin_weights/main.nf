process CREATE_THIN_WEIGHTS {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(thin_predictors_file)

    output:
    tuple val(meta), path("weights.thin"), emit: thin_weights
    tuple val("${task.process}"), val("gawk"), eval("gawk --version | sed -n '1{s/GNU Awk //;s/,.*//;p}'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    gawk '{print \$1, 1}' < ${thin_predictors_file} > weights.thin
    """

    stub:
    """
    touch weights.thin
    """
}
