process CHROMAP_INDEX {
    tag '$fasta'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::chromap=0.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        // TODO: Uncommment when singularity container becomes available
        // 'https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:ed3529ef5253d7ccbc688b6a4c5c447152685757-0' :
        'quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:ed3529ef5253d7ccbc688b6a4c5c447152685757-0' :
        'quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:ed3529ef5253d7ccbc688b6a4c5c447152685757-0' }"



    input:
    path fasta

    output:
    path "*.index"     , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = fasta.baseName
    """
    chromap \\
        -i \\
        $args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromap: \$(echo \$(chromap --version 2>&1))
    END_VERSIONS
    """
}
