process LONGRANGER_MKREF {
    tag "$meta.id"
    label 'process_medium'

    def version = '2.2.2-c1'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using longranger. Please use docker or singularity containers."
    }
    container "gitlab-registry.internal.sanger.ac.uk/tol-it/software/docker-images/longranger:${version}"

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
        longranger: \$(longranger mkref --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
