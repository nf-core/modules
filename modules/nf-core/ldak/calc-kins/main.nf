process CALC_KINS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val power
    tuple val(meta2), path(weights_file)

    output:
    tuple val(meta), path("${meta.id}.grm.bin"), path("${meta.id}.grm.id"), path("${meta.id}.grm.details"), path("${meta.id}.grm.adjust"), emit: ldak_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_prefix = "${meta.id}"
    def weights_arg = weights_file ? "--weights ${weights_file}" : ''
    def ignore_weights_arg = weights_file ? '' : "--ignore-weights YES"
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    ldak6 \\
        --calc-kins-direct ${output_prefix} \\
        --bfile ${bed.baseName} \\
        --power ${power} \\
        ${weights_arg} \\
        ${ignore_weights_arg} \\
        --max-threads ${task.cpus} ${extra_args}
    """

    stub:
    def output_prefix = "${meta.id}"
    """
    touch ${output_prefix}.grm.bin
    touch ${output_prefix}.grm.id
    touch ${output_prefix}.grm.details
    touch ${output_prefix}.grm.adjust
    """
}
