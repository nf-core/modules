process BRACKEN_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken:2.9--py38h2494328_0':
        'biocontainers/bracken:2.9--py38h2494328_0' }"

    input:
    tuple val(meta), path(kraken2db)

    output:
    tuple val(meta), path(kraken2db               , includeInputs: true), emit: db
    tuple val(meta), path("${kraken2db}/database*", includeInputs: true), emit: bracken_files
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bracken-build \\
        $args \\
        -t $task.cpus \\
        -d $kraken2db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${kraken2db}/database100mers.kmer_distrib
    touch ${kraken2db}/database100mers.kraken
    touch ${kraken2db}/database.kraken

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """
}
