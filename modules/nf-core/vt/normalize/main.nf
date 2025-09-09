process VT_NORMALIZE {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:2015.11.10--h5ef6573_4':
        'biocontainers/vt:2015.11.10--h5ef6573_4' }"

    input:
    tuple val(meta) , path(vcf), path(tbi), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("${fasta}.fai")   , emit: fai, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "-i ${intervals}" : ""

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def VERSION = "2015.11.10" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vt normalize \\
        -o ${prefix}.vcf \\
        -r ${fasta} \\
        ${regions} \\
        ${args} \\
        ${vcf}

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def VERSION = "2015.11.10" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: ${VERSION}
    END_VERSIONS
    """
}
