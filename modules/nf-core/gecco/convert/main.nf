process GECCO_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gecco:0.10.1--pyhdfd78af_0':
        'quay.io/biocontainers/gecco:0.10.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(clusters), path(gbk)
    val(mode)
    val(format)

    output:
    tuple val(meta), path("${prefix}/*.gff")        , emit: gff     , optional: true
    tuple val(meta), path("${prefix}/*.region*.gbk"), emit: bigslice, optional: true
    tuple val(meta), path("${prefix}/*.faa")        , emit: faa     , optional: true
    tuple val(meta), path("${prefix}/*.fna")        , emit: fna     , optional: true
    tuple val("${task.process}"), val('gecco'), eval("gecco -V |& sed 's/gecco //'"), emit: versions_gecco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gecco \\
        convert \\
        ${args} \\
        --jobs ${task.cpus} \\
        ${mode} \\
        --input-dir ./ \\
        --format ${format} \\
        --output ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir ${prefix}
    touch ${prefix}/${prefix}.gff
    """
}
