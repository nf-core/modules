process DRAGEN {
    tag "$meta.id"
    label 'process_long'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'registry.hub.docker.com/etycksen/dragen4:4.2.4' }" // somehow the PATH is not containing dragen, when I run it locally it is there

    input:
    tuple val(meta), path(fastq)
    path reference

    output:
    tuple val(meta), path("*.vcf"), emit: vcf, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bin_path = task.ext.bin_path ?: 'dragen'

    def input = ''
    // from fastq
    if (fastq) {
        if (fastq.size() > 2) {
            error "Error: cannot have more than 2 fastq files as input."
        } else {
            input = '-1 ' + fastq.join(' -2 ')
        }
    }

    """
    #dragen_reset

    #$bin_path \\
    echo \\
        $args \\
        $input \\
        -n $task.cpus \\
        -r $reference \\
        --output-file-prefix $prefix \\
        --output-directory \$(pwd)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """
}
