process KHMER_TRIMLOWABUND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(seq_file)

    output:
    tuple val(meta), path("${output_path}"), emit: trimmed
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_trimmed"
    if (!task.memory) {
        log.info '[KHMER_TRIMLOWABUND] Available memory not known - defaulting to 16GB. Specify process memory requirements to change this.'
        avail_mem = 16
    } else {
        avail_mem = task.memory.toGiga()
    }
    file_ext = seq_file.name - ~/\.gz$/ - ~/^[^\.]*\./
    output_path = "${prefix}.${file_ext}"
    if (args ==~ '.*--gzip.*' || args ==~ '.*--bzip.*') {
        output_path = output_path + '.gz'
    }
    if (seq_file == output_path) {
        error("Output filename is the same as input filename. Please specify a different prefix.")
    }
    """
    trim-low-abund.py \\
        $args \\
        -M ${avail_mem}e9 \\
        --output ${output_path} \\
        $seq_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( trim-low-abund.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    file_ext = seq_file.name - ~/\.gz$/ - ~/^[^\.]*\./
    output_path = "${prefix}_trimmed.${file_ext}"
    if (args ==~ '.*--gzip.*' || args ==~ '.*--bzip.*') {
        output_path = output_path + '.gz'
    }
    if (seq_file == output_path) {
        error("Output filename is the same as input filename. Please specify a different prefix.")
    }
    """
    echo "" | gzip > ${output_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( trim-low-abund.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
