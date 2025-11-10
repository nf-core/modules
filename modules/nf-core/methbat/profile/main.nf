process METHBAT_PROFILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methbat:0.16.0--h9ee0642_0':
        'biocontainers/methbat:0.16.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(files, stageAs: "inputs/*")
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.tsv"), emit: region_profile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    methbat profile \\
        --input-prefix inputs/${prefix} \\
        --input-regions ${regions} \\
        --output-region-profile ${prefix}.tsv \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | sed 's/.* //')
    END_VERSIONS
    """
}
