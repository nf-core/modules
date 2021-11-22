process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::medaka=1.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0"
    } else {
        container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"
    }

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.fa.gz"), emit: assembly
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d $assembly \\
        -o ./

    mv consensus.fasta ${prefix}.fa

    gzip -n ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
