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
        -u \\
        $args \\
        --threads $task.cpus \\

    # gzip paired reads
    if [[ -f ${reads[0]}.paired.fastq ]]; then
        gzip ${reads[0]}.paired.fastq
    fi
    if [[ -f ${reads[1]}.paired.fastq ]]; then
        gzip ${reads[1]}.paired.fastq
    fi

    # gzip unpaired reads
    if [[ -f ${reads[0]}.unpaired.fastq ]]; then
        gzip ${reads[0]}.unpaired.fastq
    fi
    if [[ -f ${reads[1]}.unpaired.fastq ]]; then
        gzip ${reads[1]}.unpaired.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
