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
    assert commands.size() > 1
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd_size  = commands.size()
    def last_args = task.ext."args$cmd_size" ?: ''
    def extension = last_args.contains("--output-fmt sam") ? "sam" :
                    last_args.contains("--output-fmt bam") ? "bam" :
                    last_args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    assert "$input" != "${prefix}.${extension}" : "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    // Compose pipe
    def cmds = commands.indexed().collect { index, cmd ->
        def first = index == 0
        def last = index == cmd_size-1
        def command = [
            "samtools $cmd",
            task.ext."args${index!=0 ? index+1 : ''}" ?: ''
        ]
        switch(cmd){
            case !"reheader":
                // The reheader has no useful option
                command << "-@ $task.cpus"
                command << fasta && last ? "--reference ${fasta}" : ''
                command << !last ? '-u' : ''
            // samtools commands have slightly different syntax
            case "collate":
                // [-o OUTPUT|-O] [INPUT|-]
                command << last ? "-o ${prefix}.${extension}" : "-O"
                command << first ? input : '-'
                break
            case ["addreplacerg", "sort", "view"]:
                // [-o OUTPUT] [INPUT|-]
                command << last ? "-o ${prefix}.${extension}" : ""
                command << first ? input : '-'
                break
            case "reheader":
                // [INPUT|-]
                command << first ? input : '-'
                break
            case ["fixmate", "markdup"]:
                // [INPUT|-] [OUTPUT|-]
                command << first ? input : '-'
                command << last ? "${prefix}.${extension}" : "-"
                break
            default:
                assert false: "$cmd is not supported"
        }
        command.join(' ')
    }

    """
    ${cmds.join(' | ')}

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
