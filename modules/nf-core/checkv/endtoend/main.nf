process CHECKV_ENDTOEND {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::checkv=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path ("${prefix}")      , emit: results
    tuple val(meta), path ("${prefix}/*.tsv"), emit: tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    checkv_db = db ? "export CHECKVDB=${db}" : ""
    """
    $checkv_db

    checkv \\
        end_to_end \\
        $args \\
        -t $task.cpus \\
        $fasta \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """
}
