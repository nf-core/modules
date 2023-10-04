process HAPLOGREP2_CLASSIFY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::haplogrep=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplogrep:2.4.0--hdfd78af_0':
        'biocontainers/haplogrep:2.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(inputfile)
    val(format)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplogrep \\
        classify \\
        $args \\
        --in $inputfile \\
        --out ${prefix}.txt \\
        --format $format

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep2: \$(echo \$(haplogrep --version 2>&1) | (sed 's/htt.*//') | (sed 's/.*v//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep2: \$(echo \$(haplogrep --version 2>&1) | (sed 's/htt.*//') | (sed 's/.*v//'))
    END_VERSIONS
    """

}
