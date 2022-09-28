process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)
    val(interleave)


    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def endedness = meta.single_end || interleave ? "-0 ${prefix}_other.fastq.gz | gzip --no-name > ${prefix}.fastq.gz" : "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -0 ${prefix}_other.fastq.gz -s ${prefix}_singleton.fastq.gz"

    // Interleaved / SE output
    // -0 ${prefix}_other.fastq.gz | gzip --no-name > ${prefix}_interleaved.fastq.gz

    // PE output
    // -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -0 ${prefix}_other.fastq.gz -s ${prefix}_singleton.fastq.gz

    // NOTE: When outputting fastq for downstream alignment, the "-n" flag should be set to avoid editing the read names
    """
    samtools \\
      collate \\
        $args \\
        --threads $task.cpus \\
        $reference \\
        -O -u \\
        $bam \\
    | samtools \\
        fastq \\
        $args \\
        --threads $task.cpus \\
        $endedness

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
