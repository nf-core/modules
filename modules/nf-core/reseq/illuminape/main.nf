process RESEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/reseq:1.1--py38hc35fec1_3':
        'biocontainers/reseq:1.1--py38hc35fec1_3' }"


    input:
    path(fasta)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: simulated_fastqgz
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reseq illuminaPE \\
        ${args} \\
        -j ${task.cpus} \\
        -r ${fasta} \\
        -b ${bam} \\
        -1 ${prefix}.1.fastq.gz \\
        -2 ${prefix}.2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reseq: \$(reseq --version |& sed 's/ReSeq version //g')
    END_VERSIONS
    """
}
