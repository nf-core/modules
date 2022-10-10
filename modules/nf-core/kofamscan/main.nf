process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kofamscan=1.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'quay.io/biocontainers/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    path profile
    path ko_list
    val format

    output:
    tuple val(meta), path('*.txt')  , optional: true, emit: txt
    tuple val(meta), path('*.tsv')  , optional: true, emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format_str = format ? "-f ${format[0]}" : ""

    switch( format ) {
        case "detail": ext = 'txt'; break
        case "detail-tsv": ext = 'tsv'; break
        case "mapper": ext = 'txt'; break
        case "mapper-oneline": ext = 'txt'; break
        default:
            ext = 'txt'
            break
    }
    """
    exec_annotation \\
        -p $profile
        -k $ko_list
        ${format_str}
        $args \\
        --cpu $task.cpus \\
        -o ${prefix}.${ext} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kofamscan: \$(echo \$(exec_annotation --version 2>&1) | sed 's/^.*exec_annotation //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
