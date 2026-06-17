process SEQKIT_TAB2FX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(text)
    val out_ext

    output:
    tuple val(meta), path("*.f*"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "fa.zst"
    """
    seqkit \\
        tab2fx \\
        ${args} \\
        --threads ${task.cpus} \\
        ${text} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "fa.zst"
    """
    echo ${args}

    touch ${prefix}.${suffix}
    """
}
