process RAVEN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::raven-assembler=1.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raven-assembler:1.6.1--h2e03b76_0' :
        'quay.io/biocontainers/raven-assembler:1.6.1--h2e03b76_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # run tool
    raven \\
        -t $task.cpus \\
        --graphical-fragment-assembly ${prefix}.gfa \\
        $args \\
        $reads | \\
        gzip -c > ${prefix}.fasta.gz

    # compress assembly graph
    gzip -c ${prefix}.gfa > ${prefix}.gfa.gz

    # get tool version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raven: \$( raven --version )
    END_VERSIONS
    """
}
