process ARTIC_GUPPYPLEX {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::artic=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.1--py_0' :
        'quay.io/biocontainers/artic:1.2.1--py_0' }"

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic \\
        guppyplex \\
        $args \\
        --directory $fastq_dir \\
        --output ${prefix}.fastq

    pigz -p $task.cpus *.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic --version 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
