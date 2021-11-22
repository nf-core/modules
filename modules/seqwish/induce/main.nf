def VERSION = '0.7.1' // Version information not provided by tool on CLI

process SEQWISH_INDUCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqwish=0.7.1' : null)

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqwish:0.7.1--h2e03b76_0' :
        'quay.io/biocontainers/seqwish:0.7.1--h2e03b76_0' }"

    input:
    tuple val(meta), path(paf), path(fasta)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    seqwish \\
        --threads $task.cpus \\
        --paf-alns=$paf \\
        --seqs=$fasta \\
        --gfa=${prefix}.gfa \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        seqwish: $VERSION
    END_VERSIONS
    """
}
