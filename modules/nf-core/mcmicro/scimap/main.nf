process MCMICRO_SCIMAP {
    tag "$meta.id"
    label 'process_single'

    container "labsyspharm/scimap:0.22.0"

    input:
    tuple val(meta), path(cellbyfeature)

    output:
    tuple val(meta), path("*.csv"), emit: annotedDataCsv, optional:true
    tuple val(meta), path("*.h5ad"), emit: annotedDataH5ad, optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    scimap-mcmicro $cellbyfeature ${prefix} ${args}
    """

}
