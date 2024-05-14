VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process MIDAS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::midas=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/midas:1.3.2--pyh5e36f6f_6':
        'biocontainers/midas:1.3.2--pyh5e36f6f_6' }"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'db/*'
    val(mode)

    output:
    tuple val(meta), path("results/*"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        run_midas.py \\
            $mode \\
            results \\
            -1 ${reads[0]} \\
            -d $db \\
            $args \\
            -t $task.cpus

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            midas: $VERSION
        END_VERSIONS
        """
    } else {
        """
        run_midas.py \\
            $mode \\
            $prefix \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -d $db \\
            $args \\
            -t $task.cpus

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            midas: $VERSION
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/species
    touch "results/species/species_profile.tsv"
    touch "results/species/log.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
    END_VERSIONS
    """
}
