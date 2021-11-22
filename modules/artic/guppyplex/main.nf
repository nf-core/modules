process ARTIC_GUPPYPLEX {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::artic=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/artic:1.2.1--py_0"
    } else {
        container "quay.io/biocontainers/artic:1.2.1--py_0"
    }

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    artic \\
        guppyplex \\
        $args \\
        --directory $fastq_dir \\
        --output ${prefix}.fastq

    pigz -p $task.cpus *.fastq
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(artic --version 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
