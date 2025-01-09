process GRIDSS_SOMATICFILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(pondir)

    output:
    tuple val(meta), path("*.high_confidence_somatic.vcf.bgz")    , emit: high_conf_sv
    tuple val(meta), path("*.all_somatic.vcf.bgz")                , emit: all_sv
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pondir_command = pondir ? "--pondir ${pondir}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gridss_somatic_filter \\
        --input $vcf \\
        ${pondir_command} \\
        --output ${prefix}.high_confidence_somatic.vcf \\
        --fulloutput ${prefix}.all_somatic.vcf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "" | bgzip > ${prefix}.high_confidence_somatic.vcf.bgz
    echo "" | bgzip > ${prefix}.all_somatic.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
