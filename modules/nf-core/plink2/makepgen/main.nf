process PLINK2_MAKEPGEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0'
        : 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.pvar"), emit: pvar
    tuple val(meta), path("*.psam"), emit: psam
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), topic: versions, emit: versions_plink2

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
