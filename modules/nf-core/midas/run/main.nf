process MIDAS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/midas:1.3.2--pyh7cba7a3_7':
        'quay.io/biocontainers/midas:1.3.2--pyh7cba7a3_7' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(db, stageAs: 'db/*')
    val(mode)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('midas'), val("1.3.2"), emit: versions_midas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        run_midas.py \\
            ${mode} \\
            ${prefix} \\
            -d ${db} \\
            ${args} \\
            -1 ${reads[0]} \\
            -t ${task.cpus}
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
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/species
    touch "${prefix}/species/species_profile.tsv"
    touch "${prefix}/species/log.txt"
    """
}
