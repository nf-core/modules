process WINDOWMASKER_CONVERT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${output}"), emit: converted
    tuple val("${task.process}"), val('windowmasker'), eval("windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker. //; s/ .*\$//'"), topic: versions, emit: versions_windowmasker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   =    task.ext.args                       ?: ''
    def prefix =    task.ext.prefix                     ?: "${meta.id}"
    def outfmt =    args.contains('-sformat binary')    ? 'binary'  :
                    args.contains('-sformat oascii')    ? 'oascii'  :
                    args.contains('-sformat obinary')   ? 'obinary' :
                    'ascii'
    output  = "${prefix}.${outfmt}"
    """
    windowmasker -convert \\
        -in ${counts} \\
        -out ${output} \\
        ${args}
    """

    stub:
    def args   =    task.ext.args                       ?: ''
    def prefix =    task.ext.prefix                     ?: "${meta.id}"
    def outfmt =    args.contains('-sformat binary')    ? 'binary'  :
                    args.contains('-sformat oascii')    ? 'oascii'  :
                    args.contains('-sformat obinary')   ? 'obinary' :
                    'ascii'
    output  = "${prefix}.${outfmt}"
    """
    touch ${output}
    """
}
