process CHROMAP_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::chromap=0.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromap:0.2.4--hd03093a_0' :
        'biocontainers/chromap:0.2.4--hd03093a_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.index"), emit: index
    path "versions.yml"              , emit: versions

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
