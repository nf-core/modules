
process BEDOPS_GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0':
        'biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    path gtf

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    cat \\
    $gtf \\
    | gtf2bed \\
    $args \\
    --attribute-key=exon_id \\
    > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtf2bed: \$(bedops --version | grep version | awk ' { print \$2 } ')
    END_VERSIONS
    """
}
