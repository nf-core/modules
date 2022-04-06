def VERSION = '0.7.2' // Version information not provided by tool on CLI

process SEQWISH_INDUCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqwish=0.7.2' : null)

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqwish:0.7.2--h2e03b76_0' :
        'quay.io/biocontainers/seqwish:0.7.2--h2e03b76_0' }"

    input:
    tuple val(meta), path(paf), path(fasta)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqwish \\
        --threads $task.cpus \\
        --paf-alns=$paf \\
        --seqs=$fasta \\
        --gfa=${prefix}.gfa \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqwish: $VERSION
    END_VERSIONS
    """
}
