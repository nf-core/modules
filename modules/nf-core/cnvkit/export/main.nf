process CNVKIT_EXPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/3935f8d7507f85fd20e073f0625bb6bdff8aa1b6044f4bbb720c9a8063ea0390/data'
:         'community.wave.seqera.io/library/cnvkit:0.9.14--288e98d6210b7304' }"

    input:
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: output
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    cnvkit.py export \\
        ${args} \\
        ${cns} \\
        -o ${prefix}.${suffix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    touch ${prefix}.${suffix}
    """
}
