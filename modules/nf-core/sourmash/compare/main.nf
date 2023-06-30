process SOURMASH_COMPARE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::sourmash=4.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.6.1--hdfd78af_0':
        'biocontainers/sourmash:4.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(signatures)
    path file_list // optional file
    val save_numpy_matrix
    val save_csv

    output:
    tuple val(meta), path("*comp.npy.labels.txt"), path("*comp.npy"), optional:true, emit: matrix
    path "*comp.csv"       , optional:true, emit: csv
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args     ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def comp   = save_numpy_matrix ? "--output comp.npy"  : ''
    def csv    = save_csv          ? "--csv comp.csv" : ''
    if ( !save_numpy_matrix && !save_csv ) error "Supply either save_numpy_matrix, save_csv, or both or no output will be created"
    def ffile = file_list ? "--from-file ${file_list}" : ''
    def sigs = signatures ? "${signatures.join(' ')}" : ''
    if ( !file_list && !signatures ) error "Supply either signatures, file_list, or both"
    """
    sourmash compare \\
        $args \\
        --processes $task.cpus \\
        ${comp} \\
        ${csv} \\
        ${ffile} \\
        ${sigs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
