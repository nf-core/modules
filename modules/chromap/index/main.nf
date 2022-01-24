process CHROMAP_INDEX {
    tag '$fasta'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::chromap=0.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromap:0.1.5--h9a82719_0' :
        'quay.io/biocontainers/chromap:0.1.5--h9a82719_0' }"


    input:
    path fasta

    output:
    path "*.index"     , emit: index
    path "versions.yml", emit: versions

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
