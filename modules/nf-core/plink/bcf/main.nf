process PLINK_BCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    tuple val(meta), path("*.bim"), emit: bim, optional: true
    tuple val(meta), path("*.fam"), emit: fam, optional: true

    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink \\
        --bcf ${bcf} \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}
    """
}
