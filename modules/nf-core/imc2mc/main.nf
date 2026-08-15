process IMC2MC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/imc2mc:0.0.3"

    input:
    tuple val(meta), path(txtfile)

    output:
    tuple val(meta), path("*.tif"), emit: tif
    tuple val("${task.process}"), val('imc2mc'), eval("imc2mc --version | sed 's/v//g'"), emit: versions_imc2mc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    imc2mc \\
        -i ${txtfile} \\
        -o "${prefix}.tif" \\
        $args

    sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' ${prefix}.tif
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tif
    """
}
