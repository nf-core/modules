process FIND_CONCATENATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pigz:2.3.4'
        : 'biocontainers/pigz:2.3.4'}"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_list = files_in.collect { it.toString() }

    // choose appropriate concatenation tool depending on input and output format

    // | input     | output     | command1 | command2 |
    // |-----------|------------|----------|----------|
    // | gzipped   | gzipped    | cat      |          |
    // | ungzipped | ungzipped  | cat      |          |
    // | gzipped   | ungzipped  | zcat     |          |
    // | ungzipped | gzipped    | cat      | pigz     |


    // get file extensions, if extension is .gz then get the second to last extension as well
    file_extensions = files_in.collect { in_file -> in_file.name - in_file.getBaseName(in_file.name.endsWith('.gz') ? 2 : 1) }.toSet()

    // Use input file ending as default
    prefix = task.ext.prefix ?: "${meta.id}${file_extensions[0]}"
    pattern_string = generatePatternString(file_extensions.toList())

    out_zip = prefix.endsWith('.gz')
    in_zip = file_list[0].endsWith('.gz')
    command1 = in_zip && !out_zip ? 'zcat' : 'cat'
    command2 = !in_zip && out_zip ? "| pigz -c -p ${task.cpus} ${args2}" : ''
    if (file_list.contains(prefix.trim())) {
        error(
            "The name of the input file can't be the same as for the output prefix in the " + "module FIND_CONCATENATE (currently `${prefix}`). Please choose a different one."
        )
    }

    """
    find . -maxdepth 1 \\( -not -name '.*'  ${pattern_string} \\) \\
         -exec sh -c "${command1} ${args} {} ${command2} >> ${prefix}" \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def file_list = files_in.collect { it.toString() }
    prefix = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.'))}"
    if (file_list.contains(prefix.trim())) {
        error(
            "The name of the input file can't be the same as for the output prefix in the " + "module FIND_CONCATENATE (currently `${prefix}`). Please choose a different one."
        )
    }
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}

def generatePatternString(fileExtensionList) {
    if (!fileExtensionList || fileExtensionList.isEmpty()) {
        return ""
    }

    if (fileExtensionList.size() == 1) {
        return "-name '*${fileExtensionList[0]}'"
    }

    def patternString = "-name '*${fileExtensionList[0]}' "
    fileExtensionList[1..-1].each {
        patternString += "-o -name '*${it}' "
    }
    return patternString.trim()
}
