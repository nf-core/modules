process HAPLOGREP3_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplogrep3:3.2.2--hdfd78af_0':
        'biocontainers/haplogrep3:3.2.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(inputfile)

    output:
    tuple val(meta)             , path("*.txt")                                                                           , emit: txt
    tuple val("${task.process}"), val('haplogrep3'), eval("haplogrep3 | sed -n 's/.*Haplogrep 3 \\([0-9.]\\+\\).*/\\1/p'"), emit: versions_haplogrep3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplogrep3 \\
        classify \\
        $args \\
        --in $inputfile \\
        --out ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """

}
