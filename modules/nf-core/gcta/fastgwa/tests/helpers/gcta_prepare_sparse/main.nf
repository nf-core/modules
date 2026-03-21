process GCTA_TEST_PREPARE_SPARSE {
    tag "${meta.id}"
    label "process_medium"
    conda "${moduleDir}/../../../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val(grm_cutoff)

    output:
    tuple val(meta), path("${meta.id}.grm.id"), path("${meta.id}.grm.sp"), emit: sparse_grm

    script:
    def bfile_prefix = bed.name.replaceFirst(/\.bed$/, "")
    def dense_prefix = "${meta.id}_dense"
    def sparse_prefix = meta.id

    """
    set -euo pipefail

    gcta \\
        --bfile ${bfile_prefix} \\
        --make-grm \\
        --out ${dense_prefix} \\
        --thread-num ${task.cpus}

    gcta \\
        --grm ${dense_prefix} \\
        --make-bK-sparse ${grm_cutoff} \\
        --out ${sparse_prefix} \\
        --thread-num ${task.cpus}
    """
}
