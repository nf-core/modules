process RACON {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::racon=1.4.20" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1"
    } else {
        container "quay.io/biocontainers/racon:1.4.20--h9a82719_1"
    }

    input:
    tuple val(meta), path(reads), path(assembly), path(paf)

    output:
    tuple val(meta), path('*_assembly_consensus.fasta.gz') , emit: improved_assembly
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    racon -t "${task.cpus}" \\
        "${reads}" \\
        "${paf}" \\
        $options.args \\
        "${assembly}" > \\
        ${prefix}_assembly_consensus.fasta

    gzip -n ${prefix}_assembly_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
