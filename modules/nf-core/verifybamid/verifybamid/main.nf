process VERIFYBAMID_VERIFYBAMID {
    tag "${meta.id}"
    label "process_single"

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid%3A1.1.3--h5b5514e_6':
        'quay.io/biocontainers/verifybamid:1.1.3--h5b5514e_6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path refvcf

    output:
    tuple val(meta), path("*.log")    , emit: log
    tuple val(meta), path("*.selfSM") , emit: selfsm , optional:true
    tuple val(meta), path("*.depthSM"), emit: depthsm, optional:true
    tuple val(meta), path("*.selfRG") , emit: selfrg , optional:true
    tuple val(meta), path("*.depthRG"), emit: depthrg, optional:true
    tuple val(meta), path("*.bestSM") , emit: bestsm , optional:true
    tuple val(meta), path("*.bestRG") , emit: bestrg , optional:true
    tuple val("${task.process}"), val('verifybamid'), eval("echo '1.1.3'"), emit: versions_verifybamid, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    verifyBamID \\
        --bam ${bam} \\
        --vcf ${refvcf} \\
        --out ${prefix} \\
        ${args} \\
        > ${prefix}.log
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.log
    touch ${prefix}.bestSM
    touch ${prefix}.bestRG
    touch ${prefix}.depthRG
    touch ${prefix}.depthSM
    touch ${prefix}.selfRG
    touch ${prefix}.selfSM
    """
}
