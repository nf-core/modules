process VT_DECOMPOSEBLOCKSUB {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d2ef4f2d3563e300c139cd56c0b6196c00c9cceca2d120384b5a0615267f630/data'
:         'community.wave.seqera.io/library/vt:0.57721--d3fc8e5bf1a4b0de' }"

    input:
    tuple val(meta), path(vcf), path(index), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    tuple val("${task.process}"), val('vt'), val('2015.11.10'), topic: versions, emit: versions_vt // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_decomposeblocksub"

    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def bed = intervals ? "-i ${intervals}" : ""
    """
    vt decompose_blocksub \\
        -o ${prefix}.vcf.gz \\
        ${bed} \\
        ${args} \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_decomposeblocksub"

    if ("$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
