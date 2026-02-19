process PDB2PQR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b84828ce9f75233b610b5778c4bb93078958515bbd8c6fe6aa77882dd22f211/data':
        'community.wave.seqera.io/library/python_pip_pdb2pqr:dfc1bd7d388ba091' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*.pqr"), emit: pqr
    tuple val(meta), path("*.in") , emit: conf
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('pdb2pqr'), eval("pdb2pqr --version | sed 's/^[^ ]* //'"), topic: versions, emit: versions_pdb2pqr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pdb2pqr \\
        $args \\
        --apbs-input=${prefix}.in \\
        ${pdb} \\
        ${prefix}.pqr
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.pqr
    touch ${prefix}.in
    touch ${prefix}.log
    """
}
