process LDAK_FILTER_RELATEDNESS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(grm_details), path(grm_adjust)

    output:
    tuple val(meta), path("${meta.id}.keep"), path("${meta.id}.lose"), path("${meta.id}.maxrel"), emit: filtered_list
    tuple val(meta), path("${meta.id}_unrel.grm.bin"), path("${meta.id}_unrel.grm.id"), path("${meta.id}_unrel.grm.details"), path("${meta.id}_unrel.grm.adjust"), emit: filtered_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filter_args = task.ext.args ?: ''
    def sub_grm_args = task.ext.args2 ?: ''
    def grm_prefix = meta.grm_prefix ?: meta.id

    """
    set -euo pipefail

    ldak6 \
        --filter ${meta.id} \
        --grm ${grm_prefix} \
        --max-threads ${task.cpus} ${filter_args}

    ldak6 \
        --sub-grm ${meta.id}_unrel \
        --grm ${grm_prefix} \
        --keep ${meta.id}.keep \
        --max-threads ${task.cpus} ${sub_grm_args}
    """

    stub:
    """
    touch ${meta.id}.keep
    touch ${meta.id}.lose
    echo "0.0" > ${meta.id}.maxrel
    touch ${meta.id}_unrel.grm.bin
    touch ${meta.id}_unrel.grm.id
    touch ${meta.id}_unrel.grm.details
    touch ${meta.id}_unrel.grm.adjust
    """
}
