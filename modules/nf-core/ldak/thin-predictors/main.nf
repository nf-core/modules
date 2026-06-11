process THIN_PREDICTORS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}.in"), emit: thin_predictors
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_prefix = "${meta.id}"
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    ldak6 \\
        --thin ${output_prefix} \\
        --bfile ${bed.baseName} \\
        --window-prune 0.98 \\
        --window-kb 100 \\
        --max-threads ${task.cpus} ${extra_args}
    """

    stub:
    def output_prefix = "${meta.id}"
    """
    touch ${output_prefix}.in
    """
}
