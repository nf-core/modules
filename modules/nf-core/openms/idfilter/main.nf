process OPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=3.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.3.0--h0656172_8' :
        'biocontainers/openms:3.3.0--h0656172_8' }"

    input:
    tuple val(meta), path(id_file), path(filter_file)

    output:
    tuple val(meta), path("*.{idXML,consensusXML}"), emit: filtered
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${id_file.getExtension()}"
    // Optional filtering via filter_file
    def filter_citerion = task.ext.args2 ?: "-whitelist:peptides"
    def filter = filter_file ? "${filter_citerion} ${filter_file}" : ""

    """
    IDFilter -in $id_file \\
        -out ${prefix}.${suffix} \\
        -threads $task.cpus \\
        $filter \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${id_file.getExtension()}"
    // Optional filtering via filter_file
    def filter_citerion = task.ext.args2 ?: "-whitelist:peptides"
    def filter = filter_file ? "${filter_citerion} ${filter_file}" : ""

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
