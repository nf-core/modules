process VERIFYBAMID_VERIFYBAMID {
    tag "${meta.id}"
    label "process_single"

    conda "bioconda::verifybamid=1.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid%3A1.1.3--h5b5514e_6':
        'biocontainers/verifybamid:1.1.3--h5b5514e_6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path refvcf

    output:
    tuple val(meta), path("*.log")                                  , emit: log
    tuple val(meta), path("*.selfSM")                , optional:true, emit: selfsm
    tuple val(meta), path("*.depthSM")               , optional:true, emit: depthsm
    tuple val(meta), path("*.selfRG")                , optional:true, emit: selfrg
    tuple val(meta), path("*.depthRG")               , optional:true, emit: depthrg
    tuple val(meta), path("*.bestSM")                , optional:true, emit: bestsm
    tuple val(meta), path("*.bestRG")                , optional:true, emit: bestrg
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    verifyBamID \\
        --bam ${bam} \\
        --vcf ${refvcf} \\
        --out ${prefix} \\
        ${args} \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid: $VERSION
    END_VERSIONS
    """
}
