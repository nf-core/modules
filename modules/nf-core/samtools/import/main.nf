process SAMTOOLS_IMPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0':
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sam") , emit: sam,  optional: true
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt bam") ? "bam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    def input = reads instanceof List && meta.single_end ? reads.join(" -0") :               // multiple single-end files
                reads instanceof List && !meta.single_end ? "-1 ${reads[0]} -2 ${reads[1]}": // paired end file
                meta.single_end ? "-0 $reads" :                                              // single single-end file
                !meta.single_end ? "-s $reads":                                              // interleave paired-end file
                reads                                                                        // if all else fails, just add the reads without flags
    """
    samtools \\
        import \\
        $input \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
