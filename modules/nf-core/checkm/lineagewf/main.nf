process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::checkm-genome=1.2.1" : null)
    def container_image = "/checkm-genome:1.2.1--pyhdfd78af_0"
                                                   container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fasta)
    val fasta_ext
    path db

    output:
    tuple val(meta), path("${prefix}")           , emit: checkm_output
    tuple val(meta), path("${prefix}/lineage.ms"), emit: marker_file
    tuple val(meta), path("${prefix}.tsv")       , emit: checkm_tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args   ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    checkm_db = db ? "export CHECKM_DATA_PATH=${db}" : ""
    """
    $checkm_db

    checkm \\
        lineage_wf \\
        -t $task.cpus \\
        -f ${prefix}.tsv \\
        --tab_table \\
        --pplacer_threads $task.cpus \\
        -x $fasta_ext \\
        $args \\
        . \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
