process SEQKIT_PAIR {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::seqkit=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.paired.fastq.gz")                  , emit: reads
    tuple val(meta), path("*.unpaired.fastq.gz"), optional: true, emit: unpaired_reads
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        pair \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        $args \\
        --threads $task.cpus

    # gzip fastq
    find . -maxdepth 1 -name "*.fastq" -exec gzip {} \;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
