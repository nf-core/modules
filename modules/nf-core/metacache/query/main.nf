process METACACHE_QUERY {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacache:2.5.0--h077b44d_0':
        'biocontainers/metacache:2.5.0--h077b44d_0' }"

    input:
    tuple val(meta), path(reads)
    path (db, stageAs: 'db/*')

    output:
    tuple val(meta), path("*.txt")          , emit: mapping_results
    tuple val(meta), path("*abundances.txt"), emit: abundances, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = meta.single_end ? reads : "${reads[0]} ${reads[1]} -pairfiles"
    """
    dbmeta=`find -L ${db}/ -name "*.meta" | head -n 1`
    metacache \\
        query \\
        \$dbmeta \\
        ${input_file} \\
        $args \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = meta.single_end ? reads : "${reads[0]} ${reads[1]} -pairfiles"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """
}
