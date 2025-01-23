// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.
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

    // Use input file ending as default
    prefix = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
    file_extensions = getSetofFileSuffixes(file_list)
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

// for .gz files also include the second to last extension if it is present. E.g., .fasta.gz
def getFileSuffix(filename) {
    def match = filename =~ /^.*?((\.\w{1,5})?(\.\w{1,5}\.gz$))/
    return match ? match[0][1] : filename.substring(filename.lastIndexOf('.'))
}

def getSetofFileSuffixes(in_files) {
    return in_files.collect { getFileSuffix(it) }.toSet()
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
