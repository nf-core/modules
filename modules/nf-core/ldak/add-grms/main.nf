process LDAK_ADD_GRMS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(mgrm_file), path(grm_files)

    output:
    tuple val(meta), path("${meta.id}.grm.bin"), path("${meta.id}.grm.id"), path("${meta.id}.grm.details"), path("${meta.id}.grm.adjust"), emit: combined_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    output_prefix = "${meta.id}"
    extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    ldak6 \
        --add-grm ${output_prefix} \
        --mgrm ${mgrm_file} \
        --max-threads ${task.cpus} ${extra_args}
    """

    stub:
    output_prefix = "${meta.id}"
    """
    touch ${output_prefix}.grm.bin
    touch ${output_prefix}.grm.id
    touch ${output_prefix}.grm.details
    touch ${output_prefix}.grm.adjust
    """
}
