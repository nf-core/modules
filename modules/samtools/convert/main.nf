process SAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(input), path(index)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.cram"), path("*.crai") , emit: cram_crai, optional: true
    tuple val(meta), path("*.bam"), path("*.bai")   , emit: bam_bai, optional:true
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        $args \\
        $input \\
        -o ${prefix}.${file_type}

    samtools index -@${task.cpus} ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
