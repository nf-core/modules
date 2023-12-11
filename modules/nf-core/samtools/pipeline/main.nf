process SAMTOOLS_PIPELINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1':
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    val commands

    output:
    tuple val(meta), path("*.{bam,cram,sam}"), emit: output
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Check that we are asked to run more than 1 command
    def n_commands = commands.size()
    if (n_commands < 2) error "SAMTOOLS_PIPELINE is used to chain 2 or more samtools commands"

    // Fetch the arguments
    def all_args = []
    for (int index = 0; index < n_commands; index++) {
        all_args.add( task.ext.args?[ commands[index] ] ?: '' )
    }
    def last_args = all_args[-1]

    // Output file
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = last_args.contains("--output-fmt sam") ? "sam" :
                    last_args.contains("--output-fmt bam") ? "bam" :
                    last_args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    def pipeline_command = ""
    for (int index = 0; index < n_commands; index++) {
        def this_command = commands[index]
        def is_first_command = (index == 0)
        def is_last_command = (index == (n_commands - 1))
        pipeline_command += """
        samtools \\
            ${this_command} \\
            ${all_args[index]} \\
            -@ $task.cpus \\
        """

        if (is_last_command) {
            // All commands support --reference, except reheader
            if (! ["reheader"].contains(this_command)) {
                pipeline_command += fasta ? " --reference ${fasta}" : ""
            }
        } else {
            // All commands support uncompressed output, except reheader
            if (! ["reheader"].contains(this_command)) {
                pipeline_command += " -u"
            }
        }
        this_input = (is_first_command ? " $input" : " -")

        // samtools commands have slightly different syntax
        if (["collate"].contains(this_command)) {
            // [-o OUTPUT|-O] [INPUT|-]
            pipeline_command += is_last_command ? " -o ${prefix}.${extension}" : " -O"
            pipeline_command += this_input

        } else if (["addreplacerg", "sort", "view"].contains(this_command)) {
            // [-o OUTPUT] [INPUT|-]
            pipeline_command += is_last_command ? " -o ${prefix}.${extension}" : ""
            pipeline_command += this_input

        } else if (["reheader"].contains(this_command)) {
            // [INPUT|-]
            pipeline_command += this_input

        } else if (["fixmate", "markdup"].contains(this_command)) {
            // [INPUT|-] [OUTPUT|-]
            pipeline_command += this_input
            pipeline_command += is_last_command ? " ${prefix}.${extension}" : " -"

        } else {
            error "${this_command} is not supported"
        }

        if (!is_last_command) {
            pipeline_command += " | \\"
        }
    }

    """
    $pipeline_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
