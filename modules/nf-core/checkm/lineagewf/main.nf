process CHECKM_LINEAGEWF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e77f70239b60110da040c4307b8048749ed1fc86262e07d27f1eb12a314d14f/data'
        : 'community.wave.seqera.io/library/checkm-genome:1.2.5--8d1d1a2477a013ce'}"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")
    val fasta_ext
    path db

    output:
    tuple val(meta), path("${prefix}"), emit: checkm_output
    tuple val(meta), path("${prefix}/lineage.ms"), emit: marker_file
    tuple val(meta), path("${prefix}.tsv"), emit: checkm_tsv
    tuple val("${task.process}"), val('checkm'), eval("checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//'"), emit: versions_checkm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def checkm_db = db ? "export CHECKM_DATA_PATH=${db}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${checkm_db}

    checkm \\
        lineage_wf \\
        -t ${task.cpus} \\
        -f ${prefix}.tsv \\
        --tab_table \\
        --pplacer_threads ${task.cpus} \\
        -x ${fasta_ext} \\
        ${args} \\
        input_bins/ \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p checkm_dummy_db
    export CHECKM_DATA_PATH=\$PWD/checkm_dummy_db
    mkdir ${prefix}/
    touch ${prefix}/lineage.ms ${prefix}.tsv
    """
}
