def VERSION="1.0.2"

process DEEPARG_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::deeparg=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeparg:1.0.2--pyhdfd78af_1' :
        'quay.io/biocontainers/deeparg:1.0.2--pyhdfd78af_1' }"
    /*
    We have to force singularity to run with -B to allow reading of a problematic file with borked read-write permissions in an upstream dependency (theanos).
        Original report: https://github.com/nf-core/funcscan/issues/23
    */
    containerOptions { "${workflow.containerEngine}" == 'singularity' ? '-B $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO' : '' }

    input:
    tuple val(meta), path(fasta), val(model)
    path(db)

    output:
    tuple val(meta), path("*.align.daa")            , emit: daa
    tuple val(meta), path("*.align.daa.tsv")        , emit: daa_tsv
    tuple val(meta), path("*.mapping.ARG")          , emit: arg
    tuple val(meta), path("*.mapping.potential.ARG"), emit: potential_arg
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    deeparg \\
        predict \\
        $args \\
        -i $fasta \\
        -o ${prefix} \\
        -d $db \\
        --model $model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: $VERSION
    END_VERSIONS
    """
}
