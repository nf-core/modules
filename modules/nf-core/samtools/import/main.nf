process SAMTOOLS_IMPORT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0':
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sam") , emit: sam,     optional: true
    tuple val(meta), path("*.bam") , emit: bam,     optional: true
    tuple val(meta), path("*.cram"), emit: cram,    optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt bam") ? "bam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    def input = reads instanceof List && meta.single_end ? reads.join(" -0") :              // multiple single-end files
                reads instanceof List && !meta.single_end ? "-1 $reads[0] -2 $reads[1]":    // paired end file
                meta.single_end ? "-0 $reads" :                                             // single single-end file
                !meta.single_end ? "-s $reads":                                             // interleave paired-end file
                reads                                                                       // if all else fails, just add the reads without flags
    """
    samtools \\
        import \\
        $input \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
