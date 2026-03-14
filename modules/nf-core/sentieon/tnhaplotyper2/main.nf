process SENTIEON_TNHAPLOTYPER2 {
    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(germline_resource)
    tuple val(meta6), path(germline_resource_tbi)
    tuple val(meta7), path(panel_of_normals)
    tuple val(meta8), path(panel_of_normals_tbi)
    val emit_orientation_data
    val emit_contamination_data

    output:
    tuple val(meta), path("*.orientation_data.tsv"),   emit: orientation_data, optional: true
    tuple val(meta), path("*.contamination_data.tsv"), emit: contamination_data, optional: true
    tuple val(meta), path("*.segments"),               emit: contamination_segments, optional: true
    tuple val(meta), path("*.stats"),                  emit: stats
    tuple val(meta), path("*.vcf.gz"),                 emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),             emit: index
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // options for "sentieon driver"
    def args2 = task.ext.args2 ?: '' // options for the TNhaplotyper2 algorithm. It could be something like "--tumor_sample <tumour_id> --normal_sample <normal_id>"
    def args3 = task.ext.args3 ?: '' // options for the OrientationBias algorithm. It could be something like "--tumor_sample <tumour_id>"
    def args4 = task.ext.args4 ?: '' // options for the ContaminationModel algorithm. It could be something like "--tumor_sample <tumour_id> --normal_sample <normal_id>"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gr_command = germline_resource ? "--germline_vcf ${germline_resource}" : ""
    def interval_command = intervals ? "--interval ${intervals}" : ""
    def pon_command = panel_of_normals ? "--pon ${panel_of_normals}" : ""
    def inputs = input.collect {in ->  "-i ${in}" }.join(" ")
    def orientation_bias_cmd = ""
    def contamination_cmd = ""

    if (emit_orientation_data) {
        orientation_bias_cmd = "--algo OrientationBias ${args3} ${prefix}.orientation_data.tsv"
    }

    if (emit_contamination_data) {
        contamination_cmd = "--algo ContaminationModel ${args4} --vcf ${germline_resource} --tumor_segments ${prefix}.segments ${prefix}.contamination_data.tsv"
    }

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        ${args} \\
        ${inputs} \\
        ${interval_command} \\
        --algo TNhaplotyper2 \\
        ${args2} \\
        ${gr_command} \\
        ${pon_command} \\
        ${prefix}.vcf.gz \\
        ${orientation_bias_cmd} \\
        ${contamination_cmd}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def orientation = emit_orientation_data ? "touch ${prefix}.orientation_data.tsv" : ""
    def contamination = emit_contamination_data ? "touch ${prefix}.{contamination_data.tsv,segments}" : ""
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    ${orientation}
    ${contamination}
    """
}
