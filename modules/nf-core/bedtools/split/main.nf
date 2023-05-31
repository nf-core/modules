process BEDTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3':
        'biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(bed), val(count)

    output:
    tuple val(meta), path("*.bed"), emit: beds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bedtools \\
        split \\
        $args \\
        -n ${count} \\
        -i ${bed} \\
        -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    create_beds = (1..count).collect {
        number = "0".multiply(4 - it.toString().size()) + "${it}"
        "    touch ${prefix}.${number}.bed"
    }.join("\n")

    """
    ${create_beds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
