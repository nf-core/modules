
process GRIDSS_GRIDSSSOMATICFILTER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'quay.io/biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bwa_index)
    tuple val(meta5), path(pondir)
    tuple val(meta6), path(rds)

    output:
    tuple val(meta), path("*.high_confidence_somatic.vcf"), emit: high_conf_sv
    tuple val(meta), path("*.all_somatic.vcf")            , emit: all_sv
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pondir = pondir ? "--pondir ${pondir}" : ""
    def opts = rds ? "--opts ${rds}" : ""
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    def ref = bwa_index ? "--ref ${fasta}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ${bwa}
    gridss_somatic_filter \\
        --input $vcf \\
        $ref \\
        $pondir \\
        $opts \\
        --output ${prefix}.high_confidence_somatic.vcf \\
        --fulloutput ${prefix}.all_somatic.vcf \\
        -n 1 \\
        -t 2 \\
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
    touch ${prefix}.high_confidence_somatic.vcf
    touch ${prefix}.all_somatic.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
