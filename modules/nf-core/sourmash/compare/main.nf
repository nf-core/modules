process SOURMASH_COMPARE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.14--hdfd78af_0':
        'biocontainers/sourmash:4.8.14--hdfd78af_0' }"

    input:
    tuple val(meta), path(signatures)
    path file_list // optional file
    val save_numpy_matrix
    val save_csv

    output:
    tuple val(meta), path("*comp.npy")            , emit: matrix, optional:true
    tuple val(meta), path("*comp.npy.labels.txt") , emit: labels, optional:true
    tuple val(meta), path("*comp.csv")            , emit: csv   , optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args     ?: ''
    def prefix = task.ext.prefix   ?: "${meta.id}"
    def comp   = save_numpy_matrix ? "--output ${prefix}_comp.npy"  : ''
    def csv    = save_csv          ? "--csv ${prefix}_comp.csv" : ''
    if ( !save_numpy_matrix && !save_csv ) error "Supply either save_numpy_matrix, save_csv, or both or no output will be created"
    def ffile = file_list ? "--from-file ${file_list}" : ''
    def sigs = signatures ? "${signatures.sort{it.toString()}.join(' ')}" : ''
    if ( !file_list && !signatures ) error "Supply either signatures, file_list, or both"
    """
    sourmash compare \\
        $args \\
        --processes ${task.cpus} \\
        ${comp} \\
        ${csv} \\
        ${ffile} \\
        ${sigs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_comp.npy.labels.txt
    touch ${prefix}_comp.npy
    touch ${prefix}_comp.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
