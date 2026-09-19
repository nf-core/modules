process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs2:2.2.9.1--py39hff71179_1':
        'quay.io/biocontainers/macs2:2.2.9.1--py39hff71179_1' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val macs2_gsize

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    tuple val(meta), path("*.gappedPeak")            , emit: gapped, optional:true
    tuple val(meta), path("*.bed")                   , emit: bed, optional:true
    tuple val(meta), path("*.bdg")                   , emit: bdg, optional:true
    tuple val("${task.process}"), val('macs2'), eval("macs2 --version | sed 's/macs2 //'"), emit: versions_macs2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    def format    = meta.single_end ? 'BAM' : 'BAMPE'
    def control   = controlbam ? "--control $controlbam" : ''
    if(args_list.contains('--format')){
        def id = args_list.findIndexOf{args_i -> args_i=='--format'}
        format = args_list[id+1]
        args_list.remove(id+1)
        args_list.remove(id)
    }
    """
    macs2 \\
        callpeak \\
        ${args_list.join(' ')} \\
        --gsize ${macs2_gsize} \\
        --format ${format} \\
        --name ${prefix} \\
        --treatment ${ipbam} \\
        ${control}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gappedPeak
    touch ${prefix}.bed
    touch ${prefix}.bdg
    touch ${prefix}.narrowPeak
    touch ${prefix}.xls
    """
}
