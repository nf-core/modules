process SEGUL_ALIGNTRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/segul:0.23.2--hc1c3326_0':
        'biocontainers/segul:0.23.2--hc1c3326_0' }"

    input:
    tuple val(meta), path(aln, stageAs: "input_dir/*")

    output:
    tuple val(meta), path("${prefix}/trimmed_alignments/*"), emit: trimmed
    tuple val(meta), path("${prefix}/trimming_summary.csv"), emit: summary
    tuple val("${task.process}"), val('segul'), eval("segul --version | sed 's/segul //'"), topic: versions, emit: versions_segul

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def data_type     = args.contains('--datatype')      ? '' : '--datatype dna'
    def input_format  = args.contains('--input-format')  ? '' : '--input-format auto'
    def output_format = args.contains('--output-format') ? '' : '--output-format fasta'
    def trim_mode     = args.contains('--missing-data') || args.contains('--pinf') ? '' : '--missing-data 0.5'
    """
    # Rename input files to .fa if they have non-standard extensions
    mkdir -p staged_input
    for f in input_dir/*; do
        name=\$(basename "\$f")
        ext="\${name##*.}"
        case "\$ext" in
            fa|fas|fasta|nex|nexus|phy|phylip) new_name="\$name" ;;
            *) new_name="\${name%.*}.fa" ;;
        esac
        ln -s "\$(realpath "\$f")" "staged_input/\$new_name"
    done

    segul align trim \\
        --dir staged_input \\
        ${input_format} \\
        ${output_format} \\
        ${data_type} \\
        --output ${prefix} \\
        ${trim_mode} \\
        ${args}
    """

    stub:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def data_type     = args.contains('--datatype')       ? '' : '--datatype dna'
    def input_format  = args.contains('--input-format')  ? '' : '--input-format auto'
    def output_format = args.contains('--output-format') ? '' : '--output-format fasta'
    def trim_mode     = args.contains('--missing-data') || args.contains('--pinf') ? '' : '--missing-data 0.5'
    """
    echo ${data_type} ${input_format} ${output_format} ${trim_mode} ${args}

    mkdir -p ${prefix}/trimmed_alignments
    touch ${prefix}/trimmed_alignments/stub.fas
    touch ${prefix}/trimming_summary.csv
    """
}
