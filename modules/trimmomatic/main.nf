process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz")         , emit: reads
    tuple val(meta), path("*.log")                   , emit: log
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "SE ${reads}" : "PE ${reads}"
    """
    trimmomatic \\
        --cores $task.cpus \\
        $args \\
        $trimmed \\
        > ${prefix}.trimmomatic.log
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic --version)
    END_VERSIONS
    """
}
