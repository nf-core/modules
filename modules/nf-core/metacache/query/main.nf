process METACACHE_QUERY {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacache:2.5.0--h077b44d_0':
        'biocontainers/metacache:2.5.0--h077b44d_0' }"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'db/*'
    val(do_abundances)

    output:
    tuple val(meta), path("*mapping.txt")   , emit: mapping_results
    tuple val(meta), path("*abundances.txt"), emit: abundances, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = meta.single_end ? reads : "${reads[0]} ${reads[1]} -pairfiles"
    def abundance_opt = do_abundances ? "-abundances ${prefix}.abundances.txt" : ''
    """
    dbmeta=`find -L db/ -name "*.meta" | head -n 1`
    [ -n "\$dbmeta" ] || { echo 'Database file "*.meta" not found!' >&2 ; exit 1 ; }
    metacache \\
        query \\
        \$dbmeta \\
        ${input_file} \\
        ${abundance_opt} \\
        ${args} \\
        -out ${prefix}.mapping.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def abundance_opt = do_abundances ? "-abundances ${prefix}.abundances.txt" : ''
    """
    touch ${prefix}.mapping.txt
    [ -n "$abundance_opt" ] && touch ${prefix}.abundances.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """
}
