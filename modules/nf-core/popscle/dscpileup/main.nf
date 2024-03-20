process POPSCLE_DSCPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0' :
        'biocontainers/popscle:0.1beta--h2c78cec_0' }"

    input:
    tuple val(meta), path(bam), path(vcf)

    output:
    tuple val(meta), path('*.cel'), emit: cel
    tuple val(meta), path('*.plp'), emit: plp
    tuple val(meta), path('*.var'), emit: var
    tuple val(meta), path('*.umi'), emit: umi
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    popscle dsc-pileup \\
        --sam $bam \\
        --vcf $vcf \\
        --out $prefix \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle dsc-pileup: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    touch ${prefix}.cel
    touch ${prefix}.var
    touch ${prefix}.plp
    touch ${prefix}.umi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle dsc-pileup: $VERSION
    END_VERSIONS
    """
}
