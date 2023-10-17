process FFQ {
    tag "${ids.size() == 1 ? ids[0] : "${ids[0]+'..'+ids[-1]}"}"
    label 'process_low'

    conda "bioconda::ffq=0.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ffq:0.2.1--pyhdfd78af_0':
        'biocontainers/ffq:0.2.1--pyhdfd78af_0' }"

    input:
    val ids

    output:
    path "*.json"      , emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def id_list = ids.sort()
    def name = id_list.size() == 1 ? ids[0] : 'metadata'
    def prefix = task.ext.prefix ?: "${name}"
    """
    ffq \\
        ${id_list.join(' ')} \\
        $args \\
        > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ffq: \$(echo \$(ffq --help 2>&1) | sed 's/^.*ffq //; s/: A command.*\$//' )
    END_VERSIONS
    """
}
