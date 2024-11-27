process SOURCEPREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourcepredict:0.5.1--pyhdfd78af_0':
        'biocontainers/sourcepredict:0.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(kraken_parse)
    path sources
    path labels
    path(taxa_sqlite, stageAs: '.etetoolkit/*')
    path(taxa_sqlite_traverse_pkl, stageAs: '.etetoolkit/*')
    val save_embedding

    output:
    tuple val(meta), path("*.embedding.sourcepredict.csv")  , optional:true, emit: embedding
    tuple val(meta), path("*.report.sourcepredict.csv")     , emit: report
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def save_embedding = save_embedding ? "-e ${prefix}.embedding.sourcepredict.csv" : ""
    """
    export NUMBA_CACHE_DIR='./tmp'
    export HOME='./'

    sourcepredict \\
        -s $sources \\
        -l $labels \\
        $args \\
        $save_embedding \\
        -t $task.cpus \\
        -o ${prefix}.report.sourcepredict.csv \\
        ${kraken_parse}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourcepredict: \$(python -c "import sourcepredict; print(sourcepredict.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.sourcepredict.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourcepredict: \$(python -c "import sourcepredict; print(sourcepredict.__version__)")
    END_VERSIONS
    """
}
