process SHIGEIFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::shigeifinder=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigeifinder:1.3.2--pyhdfd78af_0':
        'biocontainers/shigeifinder:1.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    shigeifinder \\
        $args \\
        --output ${prefix}.tsv \\
        -t $task.cpus \\
        -i $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: $VERSION
    END_VERSIONS
    """
}
