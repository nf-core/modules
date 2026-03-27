
process MACS3_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb492856efb63a7f824f0801b1386d08468cd4b7819ddc4c21e7f10e09b4fda/data':
        'community.wave.seqera.io/library/macs3:3.0.4--e0346d811b8b428e' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val   macs3_gsize

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    tuple val("${task.process}"), val('macs3'), eval("macs3 --version | sed -e 's/macs3 //'"), topic: versions, emit: versions_macs3

    tuple val(meta), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(meta), path("*.bed")       , optional:true, emit: bed
    tuple val(meta), path("*.bdg")       , optional:true, emit: bdg

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
    macs3 \\
        callpeak \\
        ${args_list.join(' ')} \\
        --gsize $macs3_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $ipbam \\
        $control
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
