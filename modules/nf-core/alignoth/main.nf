process ALIGNOTH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f879512a9cc25a7baa078e8a87e17b585464cf222f71be231e8f8a921354493/data'
        : 'community.wave.seqera.io/library/alignoth:1.4.6--049c1999033885bf'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(bed), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.html"),    emit: html,    optional: true
    tuple val(meta), path("${prefix}"),         emit: output,  optional: true
    tuple val(meta), path("${prefix}.vl.json"), emit: vl_json, optional: true
    tuple val("${task.process}"), val('alignoth'), eval("alignoth --version | sed 's/alignoth //g'"), topic: versions, emit: versions_alignoth

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = args.contains('--html')
        ? "--html > ${prefix}.html"
        : args.contains('--json') ? "> ${prefix}.vl.json" : "--output ${prefix}"
    def dir_cmd = output_cmd.contains('--output') ? "mkdir ${prefix}" : ""
    def vcf_cmd = vcf ? "--vcf ${vcf}" : ""
    def bed_cmd = bed ? "--bed ${bed}" : ""
    // --html is added with output_cmd and --json is not a flag used by the tool and is just here for easier module usage
    def args_corrected = args.replace('--html', '').replace('--json', '').trim()
    def touch_cmd = tbi ? "touch ${bai} && touch ${tbi}" : "touch ${bai}"
    """
    ${touch_cmd}
    ${dir_cmd}

    alignoth \\
        ${args_corrected} \\
        --bam-path ${bam} \\
        --reference ${fasta} \\
        ${vcf_cmd} \\
        ${bed_cmd} \\
        ${output_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = args.contains('--html')
        ? "${prefix}.html"
        : args.contains('--json') ? "${prefix}.vl.json" : ""
    def dir_cmd = args.contains('--html') || args.contains('--json') ? "" : "mkdir ${prefix}"
    """
    echo ${args}
    ${dir_cmd}
    touch ${output_cmd}
    """
}
