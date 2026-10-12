process PLINK2_PMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0'
        : 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'}"

    input:
    tuple val(meta), path(pgen, stageAs: 'input/*'), path(pvar, stageAs: 'input/*'), path(psam, stageAs: 'input/*')
    tuple val(meta2), path(pgen2, stageAs: 'merge/*'), path(pvar2, stageAs: 'merge/*'), path(psam2, stageAs: 'merge/*')

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.pvar"), emit: pvar
    tuple val(meta), path("*.psam"), emit: psam
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink2 \\
        --pgen ${pgen} \\
        --pvar ${pvar} \\
        --psam ${psam} \\
        --pmerge ${pgen2} ${pvar2} ${psam2} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        ${args} \\
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
