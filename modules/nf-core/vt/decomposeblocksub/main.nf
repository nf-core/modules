process VT_DECOMPOSEBLOCKSUB {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:2015.11.10--h5ef6573_4':
        'biocontainers/vt:2015.11.10--h5ef6573_4' }"

    input:
    tuple val(meta), path(vcf), path(index), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: \$(echo \$(vt decompose_blocksub 2>&1 | grep 'decompose_blocksub v' | sed -n 's/.*decompose_blocksub v\\([0-9]\\+\\.[0-9]\\+\\).*/\\1/p' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: \$(echo \$(vt decompose_blocksub 2>&1 | grep 'decompose_blocksub v' | sed -n 's/.*decompose_blocksub v\\([0-9]\\+\\.[0-9]\\+\\).*/\\1/p' ))
    END_VERSIONS
    """
}
