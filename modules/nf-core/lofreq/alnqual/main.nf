process LOFREQ_ALNQUAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py310h47ef89e_10' :
        'biocontainers/lofreq:2.1.5--py310h47ef89e_10' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    lofreq \\
        alnqual \\
        $args \\
        -b \\
        $bam \\
        $fasta > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
        samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
        samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
    END_VERSIONS
    """
}
