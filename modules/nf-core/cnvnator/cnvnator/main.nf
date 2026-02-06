process CNVNATOR_CNVNATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a54de0a286dcc6f73972be37104995ba96e37ef404f9d2bff294a7e8be5b7a92/data':
        'community.wave.seqera.io/library/cnvnator:0.4.1--5a467cfadbbc668d' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(root, stageAs:'input/')
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    val step // Without this parameter Nextflow can't distinguish between different steps, and resume becomes unreliable across the workflow chain

    output:
    tuple val(output_meta), path("${prefix}.root"), emit: root, optional: true
    tuple val(output_meta), path("${prefix}.tab") , emit: tab, optional: true
    tuple val("${task.process}"), val('cnvnator'), eval("cnvnator 2>&1 | sed -n '3s/CNVnator v//p'"), topic: versions, emit: versions_cnvnator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def input_cmd = bam             ? "-tree ${bam}"      : ''
    output_meta   = bam             ? meta                : meta2
    prefix        = task.ext.prefix ?: "${output_meta.id}"
    if (fasta) {
        reference = fasta.isDirectory() ? "-d ${fasta}" : "-fasta ${fasta}"
    } else {
        reference = ''
    }
    def calls_cmd = args.contains("-call") ? "> ${prefix}.tab" : ''
    cp_cmd    = root ? "cp input/* ${prefix}.root" :""
    """
    ${cp_cmd}
    cnvnator \\
        -root ${prefix}.root \\
        ${args} \\
        ${reference} \\
        ${input_cmd} \\
        ${calls_cmd}
    """

    stub:
    def args      = task.ext.args   ?: ''
    output_meta   = bam             ? meta                : meta2
    prefix        = task.ext.prefix ?: "${output_meta.id}"
    def touch_cmd = args.contains("-call") ? "touch ${prefix}.tab" : "touch ${prefix}.root"
    """
    ${touch_cmd}
    """
}
