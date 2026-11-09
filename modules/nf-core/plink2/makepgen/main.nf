process PLINK2_MAKEPGEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dfbbefb7c4dccd8c211a1c8bf94988dcd6a58a147e5c22e2553e80af907569/data'
        : 'community.wave.seqera.io/library/plink2:2.00a5.10--223badfc86e2b6e6'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.pvar"), emit: pvar
    tuple val(meta), path("*.psam"), emit: psam
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//; s/LM\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --make-pgen \\
        --memory ${mem_mb} \\
        --threads ${task.cpus} \\
        $args \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pgen
    touch ${prefix}.pvar
    touch ${prefix}.psam
    """
}
