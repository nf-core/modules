process ECHTVAR_ENCODE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87b75cb9e32b89261e8cbdca40b219a5d58fc78ebf92d5ca97c7ca23da1b9517/data':
        'community.wave.seqera.io/library/echtvar:0.2.4--e59eba33636e3aab' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(json_filters)

    output:
    tuple val(meta), path("*.zip"), emit: db
    tuple val("${task.process}"), val('echtvar'), eval("echtvar --version | sed 's/echtvar //g'"), emit: versions_echtvar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echtvar \\
        encode \\
        $args \\
        ${prefix}.zip \\
        ${json_filters} \\
        $vcf
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.zip
    """
}
