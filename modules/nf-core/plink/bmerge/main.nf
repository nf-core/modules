process PLINK_BMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1'
        : 'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1'}"

    input:
    tuple val(meta), path(bed, stageAs: 'input/*'), path(bim, stageAs: 'input/*'), path(fam, stageAs: 'input/*')
    tuple val(meta2), path(bed2, stageAs: 'merge/*'), path(bim2, stageAs: 'merge/*'), path(fam2, stageAs: 'merge/*')

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bim"), emit: bim
    tuple val(meta), path("*.fam"), emit: fam
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --bmerge ${bed2} ${bim2} ${fam2} \\
        --make-bed \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bim
    touch ${prefix}.fam
    """
}
