process BEDTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3':
        'biocontainers/bedtools:2.30.0--h468198e_3' }"

    input:
    tuple val(meta), path(input_A), path(input_B)
    path genome_file

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = genome_file ? "-g ${genome_file} -sorted" : ""
    """
    bedtools \\
        coverage \\
        $args \\
        $reference \\
        -a $input_A \\
        -b $input_B \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}
