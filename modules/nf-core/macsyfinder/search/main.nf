process MACSYFINDER_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macsyfinder:2.1.6--pyhdfd78af_0' :
        'biocontainers/macsyfinder:2.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(proteins)
    path models
    val model_names

    output:
    tuple val(meta), path("${prefix}/*")                   , emit: results
    tuple val(meta), path("${prefix}/macsyfinder.out")     , emit: stdout   , optional: true
    tuple val(meta), path("${prefix}/macsyfinder.err")     , emit: stderr   , optional: true
    tuple val(meta), path("${prefix}/all_systems.tsv")     , emit: summary  , optional: true
    tuple val(meta), path("${prefix}/all_best_solutions*") , emit: best_solutions, optional: true
    tuple val("${task.process}"), val('macsyfinder'), eval('macsyfinder --version 2>&1 | sed "1!d;s/^.*MacSyFinder //;s/ .*$//"'), topic: versions, emit: versions_macsyfinder
    tuple val("${task.process}"), val('hmmer'), eval('hmmsearch -h 2>&1 | sed "2!d;s/^# HMMER //;s/ .*$//"'), topic: versions, emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def model_arg = model_names ? "--models ${model_names}" : ""
    """
    macsyfinder \\
        --sequence-db ${proteins} \\
        --models-dir ${models} \\
        ${model_arg} \\
        --out-dir ${prefix} \\
        --worker ${task.cpus} \\
        ${args} \\
        2>| >( tee macsyfinder.err >&2 ) \\
        | tee macsyfinder.out

    mv macsyfinder.err ${prefix}/
    mv macsyfinder.out ${prefix}/

    # Remove empty error files to avoid snapshot pollution
    find ${prefix} -name "*.err" -type f -size 0 -delete
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/macsyfinder.err
    touch ${prefix}/macsyfinder.out
    touch ${prefix}/all_systems.tsv
    touch ${prefix}/all_best_solutions.txt
    """
}
