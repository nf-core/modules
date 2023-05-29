process BEDTOOLS_SHIFT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(chrom_sizes)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    bedtools \\
        shift \\
        -i $bed \\
        -g $chrom_sizes \\
        $args \\
        > ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
