def VERSION = '0.7.1'

process SEQWISH_INDUCE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqwish=0.7.1' : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqwish:0.7.1--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/seqwish:0.7.1--h2e03b76_0"
    }

    input:
    tuple val(meta), path(paf), path(fasta)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions


    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    seqwish \\
        --threads $task.cpus \\
        --paf-alns=$paf \\
        --seqs=$fasta \\
        --gfa=${prefix}.gfa \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
