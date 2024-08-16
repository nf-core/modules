process CAFE {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cafe:5.1.0--h43eeafb_0':
        'biocontainers/cafe:5.1.0--h43eeafb_0' }"

    input:
    path(infile)
    path(tree)

    output:
    path("Out_folder")    , emit: cafe
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cafe5
        -i ${infile} \\
        -t ${tree} \\
        $args \\
        --cores ${task.cpus} \\
        -o Out_folder


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '5.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cafe: $VERSION
    END_VERSIONS
    """
}
