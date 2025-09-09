process MIDAS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/midas:1.3.2--pyh7cba7a3_7':
        'biocontainers/midas:1.3.2--pyh7cba7a3_7' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(db, stageAs: 'db/*')
    val(mode)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    if (meta.single_end) {
        """
        run_midas.py \\
            ${mode} \\
            ${prefix} \\
            -d ${db} \\
            ${args} \\
            -1 ${reads[0]} \\
            -t ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            midas: ${VERSION}
        END_VERSIONS
        """
    } else {
        """
        run_midas.py \\
            ${mode} \\
            ${prefix} \\
            -d ${db} \\
            ${args} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            midas: ${VERSION}
        END_VERSIONS
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p ${prefix}/species
    touch "${prefix}/species/species_profile.tsv"
    touch "${prefix}/species/log.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: ${VERSION}
    END_VERSIONS
    """
}
