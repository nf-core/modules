process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::checkm-genome=1.1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/checkm-genome:1.1.3--py_1"
    } else {
        container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
    }

    input:
    tuple val(meta), path(fasta)
    val fasta_ext

    output:
    tuple val(meta), path("${prefix}")    , emit: checkm_output
    tuple val(meta), path("${prefix}.tsv"), emit: checkm_tsv
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
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
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
