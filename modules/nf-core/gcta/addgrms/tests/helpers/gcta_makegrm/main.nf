process GCTA_TEST_MAKEGRM {
    tag "${meta.id}"
    label "process_medium"
    conda "${moduleDir}/../../../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(extract_variants)

    output:
    tuple val(meta), path("${meta.id}.part_1_1.grm.id"), path("${meta.id}.part_1_1.grm.bin"), path("${meta.id}.part_1_1.grm.N.bin"), emit: grm_files

    script:
    def bfile_prefix = bed.name.replaceFirst(/\.bed$/, "")
    def out_prefix = "${meta.id}.part_1_1"

    """
    set -euo pipefail

    gcta \\
        --bfile ${bfile_prefix} \\
        --extract ${extract_variants} \\
        --make-grm \\
        --out ${out_prefix} \\
        --thread-num ${task.cpus}
    """

    stub:
    def out_prefix = "${meta.id}.part_1_1"

    """
    touch ${out_prefix}.grm.id
    touch ${out_prefix}.grm.bin
    touch ${out_prefix}.grm.N.bin
    """
}
