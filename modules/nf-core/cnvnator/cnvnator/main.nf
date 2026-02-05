process CNVNATOR_CNVNATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a54de0a286dcc6f73972be37104995ba96e37ef404f9d2bff294a7e8be5b7a92/data':
        'community.wave.seqera.io/library/cnvnator:0.4.1--5a467cfadbbc668d' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(root)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(output_meta), path("${prefix}${step}.root"), emit: root, optional: true
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
    step = args.contains("-call") ? "_call"
            : args.contains("-stat") ? "_stat"
            : args.contains("-his") ? "_his"
            : args.contains("-partition") ? "_partition" : "_rd"
    def calls_cmd = args.contains("-call") ? "> ${prefix}.tab" : ''
    def cp_cmd    = root ? "cp ${root} input.root" :""
    def mv_cmd    = args.contains("-call") ? "" : "mv input.root ${prefix}${step}.root"
    """
    ${cp_cmd}
    cnvnator \\
        -root input.root \\
        ${args} \\
        ${reference} \\
        ${input_cmd} \\
        ${calls_cmd}
    ${mv_cmd}
    """

    stub:
    def args      = task.ext.args   ?: ''
    output_meta   = bam             ? meta                : meta2
    prefix        = task.ext.prefix ?: "${output_meta.id}"
    step = args.contains("-call") ? "_call"
            : args.contains("-stat") ? "_stat"
            : args.contains("-his") ? "_his"
            : args.contains("-partition") ? "_partition" : "_rd"
    def touch_cmd = args.contains("-call") ? "touch ${prefix}.tab" : "touch ${prefix}${step}.root"
    """
    ${touch_cmd}
    """
}
