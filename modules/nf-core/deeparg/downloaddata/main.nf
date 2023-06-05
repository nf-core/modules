process DEEPARG_DOWNLOADDATA {
    label 'process_single'

    conda "bioconda::deeparg=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeparg:1.0.2--pyhdfd78af_1' :
        'biocontainers/deeparg:1.0.2--pyhdfd78af_1' }"
    /*
    We have to force singularity to run with -B to allow reading of a problematic file with borked read-write permissions in an upstream dependency (theanos).
    Original report: https://github.com/nf-core/funcscan/issues/23
    */
    containerOptions { "${workflow.containerEngine}" == 'singularity' ? '-B $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO' : '' }


    input:

    output:
    path "db/"          , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION='1.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    deeparg \\
        download_data \\
        $args \\
        -o db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: $VERSION
    END_VERSIONS
    """
}
