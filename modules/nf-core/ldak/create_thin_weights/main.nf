process CREATE_THIN_WEIGHTS {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
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
