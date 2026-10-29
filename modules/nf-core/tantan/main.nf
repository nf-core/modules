process TANTAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1a8043d5f5b80b4b7af4c60f5aa09303ab22d76807542d12b4da5e5ede48678/data'
:         'community.wave.seqera.io/library/tantan:51--ad7ba29d3afed598' }"

    input:
    tuple val(meta), path(fastx)
    val(output_format)


    output:
    tuple val(meta), path("*.fasta.gz"), emit: masked_fasta,  optional: true
    tuple val(meta), path("*.txt"), emit: probabilities, optional: true
    tuple val(meta), path("*.counts.txt"), emit: counts, optional: true
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    tuple val(meta), path("*.repeats.txt"), emit: tandem_repeats, optional: true
    tuple val("${task.process}"), val('tantan'), eval("tantan --version | sed -e 's/tantan //g'"), topic: versions, emit: versions_tantan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.id
    def format = output_format?.toString()

    if (!(format in ['0', '1', '2', '3', '4'])) {
        error(
            "Invalid tantan output format '${output_format}'. " +
            "Expected one of: 0, 1, 2, 3, or 4."
        )
    }

    /*
     * output_format is provided as a process input.
     * Do not also specify -f in task.ext.args.
     */

    def format_matches = args.findAll(/(?<!\S)-f\s*([0-4])(?!\S)/)

    if (format_matches) {
        error(
            "Do not specify '-f' in task.ext.args when output_format " +
            "is provided as a process input."
        )
    }


    def output_files = [
        '0': "${prefix}.fasta",
        '1': "${prefix}.txt",
        '2': "${prefix}.counts.txt",
        '3': "${prefix}.bed",
        '4': "${prefix}.repeats.txt"
    ]

    def output_file = output_files[format]

    def gzip_command = format == '0' ? """
        gzip -c ${output_file} > ${output_file}.gz
        rm ${output_file}
    """ : ''


    """
    zcat ${fastx} |
    tantan \\
        -f ${format} \\
        $args \\
        > ${output_file}

    $gzip_command
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def format = output_format?.toString()

    if (!(format in ['0', '1', '2', '3', '4'])) {
        error(
            "Invalid tantan output format '${output_format}'. " +
            "Expected one of: 0, 1, 2, 3, or 4."
        )
    }

    def format_matches = args.findAll(/(?<!\S)-f\s*[0-4](?!\S)/)

    if (format_matches) {
        error(
            "Do not specify '-f' in task.ext.args when output_format " +
            "is provided as a process input."
        )
    }

    def output_files = [
        '0': "${prefix}.fasta.gz",
        '1': "${prefix}.txt",
        '2': "${prefix}.counts.txt",
        '3': "${prefix}.bed",
        '4': "${prefix}.repeats.txt"
    ]

    def final_output = output_files[format]

    if (final_output == fastx.toString()) {
        error(
            'Input and output names are the same; ' +
            'use "task.ext.prefix" to disambiguate.'
        )
    }

    def create_command = final_output.endsWith('.gz')
        ? "echo '' | gzip > ${final_output}"
        : "touch ${final_output}"

    """
    echo ${args}
    echo ${args2}

    ${create_command}
    """
}
