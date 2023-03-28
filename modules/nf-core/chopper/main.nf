process CHOPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::chopper=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.3.0--hd03093a_0':
        'quay.io/biocontainers/chopper:0.3.0--hd03093a_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta)                 , emit: meta
    path "versions.yml"             , emit: versions
    tuple path("*.filtered.fastq.gz")  , emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $fastq | \\
    chopper \\
        --threads $task.cpus \\
        $args | \\
    gzip >  ${prefix}.filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
