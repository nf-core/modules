process SEQSERO2 {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqsero2=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqsero2:1.2.1--py_0"
    } else {
        container "quay.io/biocontainers/seqsero2:1.2.1--py_0"
    }

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("results/*_log.txt")   , emit: log
    tuple val(meta), path("results/*_result.tsv"), emit: tsv
    tuple val(meta), path("results/*_result.txt"), emit: txt
    path "versions.yml"                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    SeqSero2_package.py \\
        $args \\
        -d results/ \\
        -n $prefix \\
        -p $task.cpus \\
        -i $seqs

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}
