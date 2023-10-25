process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda 'modules/nf-core/cat/cat/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_list = files_in.collect { it.toString() }

    // construct output file name using meta.id and file extension of the first input file
    if ( file_list[0] =~ /^.+?(\.[\w\d]{1,5}\.gz$)/) {
        // for .gz files also include the second to last extension if it is present. E.g., .fasta.gz
        prefix   = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.', file_list[0].lastIndexOf('.')-1))}"
    } else {
        // for regular file extensions, only use a single suffix, e.g. .txt, .gz, etc.
        prefix   = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.'))}"
    }

    // choose appropriate concatenation tool depending on input and output format

    // | input     | output     | command1 | command2 |
    // |-----------|------------|----------|----------|
    // | gzipped   | gzipped    | cat      |          |
    // | ungzipped | ungzipped  | cat      |          |
    // | gzipped   | ungzipped  | zcat     |          |
    // | ungzipped | gzipped    | cat      | pigz     |

    out_zip  = prefix.endsWith('.gz')
    in_zip   = file_list[0].endsWith('.gz')
    command1 = (in_zip && !out_zip) ? 'zcat' : 'cat'
    command2 = (!in_zip && out_zip) ? "| pigz -c -p $task.cpus $args2" : ''
    """
    $command1 \\
        $args \\
        ${file_list.join(' ')} \\
        $command2 \\
        > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def file_list = files_in.collect { it.toString() }
    prefix   = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.'))}"
    """
    touch $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
