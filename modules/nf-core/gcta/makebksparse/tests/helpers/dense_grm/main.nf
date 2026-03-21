process GCTA_TEST_DENSE_GRM {
    tag "${meta.id}"
    label "process_medium"
    conda "${projectDir}/modules/nf-core/gcta/makebksparse/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}.grm.id"), path("${meta.id}.grm.bin"), path("${meta.id}.grm.N.bin"), emit: dense_grm

    script:
    def bfile_prefix = bed.baseName
    """
    set -euo pipefail

    gcta \\
        --bfile "${bfile_prefix}" \\
        --make-grm \\
        --out "${meta.id}" \\
        --thread-num ${task.cpus}
    """

    stub:
    """
    touch "${meta.id}.grm.id"
    touch "${meta.id}.grm.bin"
    touch "${meta.id}.grm.N.bin"
    """
}
