process PRINSEQPLUSPLUS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::prinseq-plus-plus=1.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prinseq-plus-plus:1.2.3--hc90279e_1':
        'quay.io/biocontainers/prinseq-plus-plus:1.2.3--hc90279e_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_good_out*.fastq.gz")                  , emit: good_reads
    tuple val(meta), path("*_single_out*.fastq.gz"), optional: true, emit: single_reads
    tuple val(meta), path("*_bad_out*.fastq.gz")   , optional: true, emit: bad_reads
    tuple val(meta), path("*.log")                                 , emit: log
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        prinseq++ \\
            -threads $task.cpus \\
            -fastq ${reads} \\
            -out_name ${prefix} \\
            -out_gz \\
            -VERBOSE 1 \\
            $args \\
            | tee ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            prinseqplusplus: \$(echo \$(prinseq++ --version | cut -f 2 -d ' ' ))
        END_VERSIONS
        """
    } else {
        """
        prinseq++ \\
            -threads $task.cpus \\
            -fastq ${reads[0]} \\
            -fastq2 ${reads[1]} \\
            -out_name ${prefix} \\
            -out_gz \\
            -VERBOSE 1 \\
            $args \\
            | tee ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            prinseqplusplus: \$(echo \$(prinseq++ --version | cut -f 2 -d ' ' ))
        END_VERSIONS
        """
    }
}
