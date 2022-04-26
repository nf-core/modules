process LONGRANGER_MKREF {
    tag "$meta.id"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using longranger"
    }
    if ( workflow.containerEngine == 'singularity' || \
            workflow.containerEngine == 'docker' ) {
        exit 1, "Longranger can not be run in container environment"
    }

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("refdata-*"),           emit: folder
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    longranger mkref $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    longranger: \$(echo \$(longranger mkref --version) | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
