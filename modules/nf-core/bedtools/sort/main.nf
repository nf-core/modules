process BEDTOOLS_SORT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(intervals)
    path genome_file

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def genome_cmd = genome_file     ?  "-g $genome_file" : ""
    extension      = task.ext.suffix ?: intervals.extension
    if ("$intervals" == "${prefix}.${extension}") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    """
    bedtools \\
        sort \\
        -i $intervals \\
        $genome_cmd \\
        $args \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
