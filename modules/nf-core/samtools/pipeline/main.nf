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
    def reference = fasta ? "--reference ${fasta}" : ""
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
        if (this_command != "fixmate") {
            pipeline_command += "    -T ${prefix}.tmp.${this_command} \\\n"
        }
        pipeline_command += "    "
        if (!is_last_command) {
            pipeline_command += "-u "
        }
        this_input = (is_first_command ? "$input" : "-")
        if (["collate", "sort"].contains(this_command)) {
            if (is_last_command) {
                pipeline_command += "-o ${prefix}.${extension} "
            } else {
                if (this_command != "sort") {
                    pipeline_command += "-O "
                }
            }
            pipeline_command += this_input
        } else {
            // e.g. fixmate, markdup
            pipeline_command += this_input
            if (is_last_command) {
                pipeline_command += " ${prefix}.${extension}"
            } else {
                pipeline_command += " -"
            }
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
