process HAPLOGREP3_CLASSIFY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplogrep3:3.2.2--hdfd78af_0':
        'biocontainers/haplogrep3:3.2.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(inputfile)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: \$(echo \$(haplogrep3 2>&1) | (sed '2!d') | (sed 's/Haplogrep 3 //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: \$(echo \$(haplogrep3 2>&1) | (sed '2!d') | (sed 's/Haplogrep 3 //'))
    END_VERSIONS
    """

}
