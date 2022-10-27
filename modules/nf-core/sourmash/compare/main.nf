process SOURMASH_COMPARE {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::sourmash=4.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.5.0--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    path signatures
    val save_numpy_matrix
    val save_csv

    output:
    path "*comp"           , optional:true, emit: numpy_matrix
    path "*comp.labels.txt", optional:true, emit: labels
    path "*comp.csv"       , optional:true, emit: csv
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args     ?: ''
    def comp   = save_numpy_matrix ? "--output comp"  : ''
    def csv    = save_csv          ? "--csv comp.csv" : ''
    """
    sourmash compare \\
        $args \\
        --processes $task.cpus \\
        ${comp} \\
        ${csv} \\
        ${signatures.join(' ')}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
