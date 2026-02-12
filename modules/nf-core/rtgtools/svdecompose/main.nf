process RTGTOOLS_SVDECOMPOSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(input), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    tuple val("${task.process}"), val('rgtools'), eval("rtg version | head -n 1 | sed 's/Product: RTG Tools //'"), topic: versions, emit: versions_rtgtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ""
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def index_vcf = tbi ? "" : "rtg index ${input}"
    def avail_mem = task.memory.toGiga() + "G"
    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ${index_vcf}

    rtg RTG_MEM=$avail_mem svdecompose \\
        ${args} \\
        --input=${input} \\
        --output=${prefix}.vcf.gz

    rtg index ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo | gzip -n > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
