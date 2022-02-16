process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::macs2=2.2.7.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py38h4a8c8d9_3' :
        'quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val   macs2_gsize

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    path  "versions.yml"                             , emit: versions

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
        def id = args_list.findIndexOf{it=='--format'}
        format = args_list[id+1]
        args_list.remove(id+1)
        args_list.remove(id)
    }
    """
    macs2 \\
        callpeak \\
        ${args_list.join(' ')} \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $ipbam \\
        $control

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
    END_VERSIONS
    """
}
