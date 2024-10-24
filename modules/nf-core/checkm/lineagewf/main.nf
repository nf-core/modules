process CHECKM_LINEAGEWF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.2.3--pyhdfd78af_0' :
        'biocontainers/checkm-genome:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")
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
        input_bins/ \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
