process FIND_UNPIGZ {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pigz:2.8'
        : 'biocontainers/pigz:2.8'}"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}.*"), emit: file_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def args3 = task.ext.args3 ?: ""
    def args4 = task.ext.args4 ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    file_extensions = files_in.collect { in_file -> in_file.name - in_file.getBaseName(in_file.name.endsWith('.gz') ? 2 : 1) }.toSet()

    file_names = files_in.collect { it.toString() }

    pattern_string = generatePatternString(file_extensions.toList())

    if (!file_extensions.every { it.endsWith(".gz") }) {
        error("All files provided to this module must be gzipped (and have the .gz extension).")
    }

    if (file_names.any { it.startsWith("${prefix}") }) {
        error("No input files can start with the same name as the output prefix in the module FIND_UNPIGZ (currently '${prefix}'). Please choose a different one.")
    }

    """
    find . -maxdepth 1 \\( -not -name '.*'  ${pattern_string} \\) ${args} |\\
        sed ${args2} 's:^./::g' | sed ${args3} 's/.gz\$//g' | xargs -I{} sh -c "unpigz -cd --processes ${task.cpus} ${args4} {}.gz > ${prefix}.{}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find: \$( find --version | head -n 1 | sed 's/find (GNU findutils) //g' )
        sed: \$( sed --version | head -n 1 | sed 's/sed (GNU sed) //g' )
        xargs: \$( xargs --version | head -n 1 | sed 's/xargs (GNU findutils) //g' )
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${files_in[0].dropRight(3)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find: \$( find --version | head -n 1 | sed 's/find (GNU findutils) //g' )
        sed: \$( sed --version | head -n 1 | sed 's/sed (GNU sed) //g' )
        xargs: \$( xargs --version | head -n 1 | sed 's/xargs (GNU findutils) //g' )
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
