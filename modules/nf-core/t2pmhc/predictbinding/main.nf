process T2PMHC_PREDICTBINDING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6d701534f0619274343a2eb556136d7a61d34cae4cc89ae9c15aa34337da3b2/data'
:         'community.wave.seqera.io/library/pip_python_t2pmhc:3066a62c82b8945d' }"

    input:
    tuple val(meta), path(samplesheet), path(graphs)
    val mode

    output:
    tuple val(meta), path("*_predicted.tsv"), emit: pred
    tuple val("${task.process}"), val('t2pmhc'), eval("pip show t2pmhc | grep -i '^Version:' | sed 's/^Version: //'"), topic: versions, emit: versions_t2pmhc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    t2pmhc t2pmhc-predict-binding \\
        --mode $mode \\
        --samplesheet ${samplesheet} \\
        --saved_graphs ${graphs} \\
        --out ${prefix}_${mode}_predicted.tsv \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_${mode}_predicted.tsv
    """
}
