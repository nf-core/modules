process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::checkm-genome=1.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.1.3--py_1' :
        'quay.io/biocontainers/checkm-genome:1.1.3--py_1' }"

    input:
    tuple val(meta), path(fasta)
    val fasta_ext

    output:
    tuple val(meta), path("${prefix}")    , emit: checkm_output
    tuple val(meta), path("${prefix}.tsv"), emit: checkm_tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
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
