process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cooltools=0.5.4 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:4a5b143bd2e45d383495ba1cd27ab81d85a5440b-0' :
        'biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:4a5b143bd2e45d383495ba1cd27ab81d85a5440b-0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*tsv"), emit:tsv
    tuple val(meta), path("*.bw"), optional: true, emit: bigwig
    path("versions.yml"), emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools insulation \\
        -p ${task.cpus} \\
        -o ${prefix}_insulation.tsv \\
        ${cool} \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | grep version | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
