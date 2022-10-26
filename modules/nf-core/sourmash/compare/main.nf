process SOURMASH_COMPARE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::sourmash=4.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.5.0--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(signatures)
    val save_numpy_matrix
    val save_csv

    output:
    tuple val(meta), path("*.comp")            , optional:true, emit: numpy_matrix
    tuple val(meta), path("*.comp.labels.txt") , optional:true, emit: labels
    tuple val(meta), path("*.comp.csv")        , optional:true, emit: csv
    path "versions.yml"                        , optional:true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args     ?: ''
    def prefix = task.ext.prefix   ?: "${meta.id}"
    def comp   = save_numpy_matrix ?: "--output ${prefix}.comp"  : '' 
    def csv    = save_csv          ?: "--csv ${prefix}.comp.csv" : ''
    
    """
    sourmash compare \\
        $args \\
        --processes $task.cpus \\
        ${comp} \\
        ${csv} \\
        ${signatures}
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' ) 
    END_VERSIONS
    """
}
