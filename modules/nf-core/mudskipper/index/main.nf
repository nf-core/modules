process MUDSKIPPER_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mudskipper:0.1.0--h9f5acd7_1':
        'quay.io/biocontainers/mudskipper:0.1.0--h9f5acd7_1' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}/"), emit: index
    tuple val("${task.process}"), val('mudskipper'), eval("mudskipper -V 2>&1 | sed 's/.*mudskipper //'"), emit: versions_mudskipper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export RUST_BACKTRACE=full
    mudskipper \\
        index \\
        --gtf ${gtf} \\
        --dir-index ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}/
    touch ${prefix}/gtf.{exon,len,map,name}
    """
}
