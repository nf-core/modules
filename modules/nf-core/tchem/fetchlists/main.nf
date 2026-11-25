process TCHEM_FETCHLISTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c61b3fd2a2c757ebba98edeec47abf69e8520970fbc5fbda0a17ceb99e59ec2/data'
        : 'community.wave.seqera.io/library/tussl-chemistry:0.1.0--6043b7490da0bfc8'}"

    input:
    path cellranger // optional: a Cell Ranger install directory, for the lists only it ships (Flex); [] when none

    output:
    path "whitelists", emit: whitelists
    tuple val("${task.process}"), val('tchem'), eval("tchem --version | sed 's/^tchem //'"), emit: versions_tchem, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The barcode lists belong to their vendors and are never shipped: the script downloads each
    // one and refuses any whose barcode count, length or sha256 differs from the registry's.
    def args = task.ext.args ?: ''
    def from_cellranger = cellranger ? "--cellranger ${cellranger}" : ''
    """
    tchem-fetch-whitelists whitelists ${from_cellranger} ${args}
    tchem --whitelists whitelists --pack
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    mkdir whitelists
    touch whitelists/10x_v3_permit.txt
    """
}
