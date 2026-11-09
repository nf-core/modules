process PLINK_HOMOZYG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h031d066_5' :
        'quay.io/biocontainers/plink:1.90b6.21--h031d066_5' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.hom")         , emit: hom
    tuple val(meta), path("*.hom.indiv")   , emit: hom_indiv
    tuple val(meta), path("*.hom.summary") , emit: hom_summary
    tuple val(meta), path("*.log")         , emit: log
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink \\
        --bfile ${bed.baseName} \\
        --homozyg \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}_ROH
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ROH.hom ${prefix}_ROH.hom.indiv ${prefix}_ROH.hom.summary ${prefix}_ROH.log
    """
}
