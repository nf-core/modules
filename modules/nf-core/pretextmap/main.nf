
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61%3A44321ab4d64f0b6d0c93abbd1406369d1b3da684-0':
        'biocontainers/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:44321ab4d64f0b6d0c93abbd1406369d1b3da684-0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.pretext")  , emit: pretext
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: '' // PretextMap args
    def args2       = task.ext.args2    ?: '' // Samtools view args
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def reference   = fasta             ? "--reference ${fasta}" : ""
    """
    if [[ $input == *.pairs.gz ]]; then
        zcat $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    else
        samtools \\
            view \\
            $args2 \\
            $reference \\
            -h \\
            $input | \\
        PretextMap \\
            $args \\
            -o ${prefix}.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | sed '/Version/!d; s/.*Version //')
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | sed '/Version/!d; s/.*Version //')
        samtools: \$(samtools --version | sed '1!d; s/samtools //')
    END_VERSIONS
    """
}
