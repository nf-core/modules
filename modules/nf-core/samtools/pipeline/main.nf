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

    // Check that we are asked to run a pipeline
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
    for (int index = 0; index < commands.size(); index++) {
        def this_command = """
        samtools \\
            ${commands[index]} \\
            ${all_args[index]} \\
            -@ $task.cpus \\
        """
        if (commands[index] != "fixmate") {
            this_command += "    -T ${prefix}.tmp.${commands[index]} \\\n"
        }
        this_command += "    "
        if (index < (n_commands - 1)) {
            this_command += "-u "
        }
        this_input = (index ? "-" : "$input")
        if (["collate", "sort"].contains(commands[index])) {
            if (index == (n_commands - 1)) {
                this_command += "-o ${prefix}.${extension} "
            } else {
                if (commands[index] != "sort") {
                    this_command += "-O "
                }
            }
            this_command += this_input
        } else {
            // e.g. fixmate, markdup
            this_command += this_input
            if (index == (n_commands - 1)) {
                this_command += " ${prefix}.${extension}"
            } else {
                this_command += " -"
            }
        }
        pipeline_command += this_command
        if (index < (n_commands - 1)) {
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
}
