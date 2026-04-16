process PHENOIMAGER2MC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/phenoimager2mc:v0.3.2"

    input:
    tuple val(meta), path(tiles, stageAs: "tiles/*")

    output:
    tuple val(meta), path("*.tif"), emit: tif
    tuple val("${task.process}"), val('phenoimager2mc'), eval("phenoimager2mc --version"), topic: versions, emit: versions_phenoimager2mc
    tuple val("${task.process}"), val('python'), eval('python --version | sed -e "s/Python //g"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('ome_types'), eval('python -m pip show ome_types | sed -n "s/Version: //p"'), topic: versions, emit: versions_ome_types

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    phenoimager2mc \
        -i ${tiles} \
        -o "${prefix}.tif" \
        ${args}

    sed -i -E \
        -e 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' \
        ${prefix}.tif
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.tif
    """
}
