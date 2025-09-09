process DEEPARG_DOWNLOADDATA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deeparg:1.0.4--pyhdfd78af_0'
        : 'biocontainers/deeparg:1.0.4--pyhdfd78af_0'}"

    /*
    We have to force docker/singularity to mount a fake file to allow reading of a problematic file with borked read-write permissions in an upstream dependency (theanos).
    Original report: https://github.com/nf-core/funcscan/issues/23
    */
    containerOptions {
        ['singularity', 'apptainer'].contains(workflow.containerEngine)
            ? '-B $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO'
            : "${workflow.containerEngine}" == 'docker'
                ? '-v $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO'
                : ''
    }

    output:
    path "db/", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.0.4'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    # Theano needs a writable space and uses the home directory by default,
    # but the latter is not always writable, for instance when Singularity
    # is run in --no-home mode
    mkdir -p theano
    export THEANO_FLAGS="base_compiledir=\$PWD/theano"

    deeparg \\
        download_data \\
        ${args} \\
        -o db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '1.0.4'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "deeparg \\
        download_data \\
        ${args} \\
        -o db/"

    mkdir db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: ${VERSION}
    END_VERSIONS
    """
}
