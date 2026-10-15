process PLINK2_VCF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/278e4e95ad25ec7048ab124f16a7cf8b8141ced9e05c51f61803567cd9017f4d/data'
:         'community.wave.seqera.io/library/plink2:2.00a2.3--1c05b6e0948ac0c2' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.pgen")    , emit: pgen
    tuple val(meta), path("*.psam")    , emit: psam
    tuple val(meta), path("*.pvar")    , emit: pvar
    tuple val(meta), path("*.pvar.zst"), emit: pvar_zst, optional: true
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        $args \\
        --vcf $vcf \\
        --out ${prefix}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar
    touch ${prefix}.pvar.zst
    """
}
