process OPENMS_PEAKPICKERHIRES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.1.0--h8964181_3' :
        'biocontainers/openms:3.1.0--h8964181_3' }"

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    PeakPickerHiRes \\
        -in $mzml \\
        -out ${prefix}.mzML \\
        -threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.mzML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
