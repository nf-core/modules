process NANOMONSV_GET {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nanomonsv:0.9.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/nanomonsv:0.9.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(tumor_parse_files, arity: '8')
    tuple val(meta2), path(control_bam), path(control_bai), path(control_parse_files, arity: '0..8')
    tuple val(meta3), path(ref), path(ref_index)
    path simple_repeat_bed
    path simple_repeat_bed_index
    path control_panel_files, arity: '0..8'

    output:
    tuple val(meta), path("${prefix}.nanomonsv.result.txt"), emit: result_txt
    tuple val(meta), path("${prefix}.nanomonsv.result.vcf"), emit: result_vcf
    tuple val(meta), path("${prefix}.nanomonsv.supporting_read.txt"), emit: supporting_read
    tuple val(meta), path("${prefix}.nanomonsv.sbnd.result.txt"), emit: sbnd_result_txt, optional: true
    tuple val(meta), path("${prefix}.nanomonsv.sbnd.result.vcf"), emit: sbnd_result_vcf, optional: true
    tuple val("${task.process}"), val('nanomonsv'), eval("nanomonsv --version 2>&1 | sed 's/^nanomonsv //'"), topic: versions, emit: versions_nanomonsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix  = task.ext.prefix  ?: "${meta.id}"
    prefix2 = task.ext.prefix2 ?: "${meta2.id}"

    def simple_repeat_arg = simple_repeat_bed ? "--simple_repeat_bed ${simple_repeat_bed[0]}" : ""
    def control_bam_arg = control_bam ? "--control_bam ${control_bam[0]}" : ""

    """
    nanomonsv get \\
        ${args} \\
        --processes ${task.cpus} \\
        ${control_bam_arg} \\
        --control_prefix ${prefix2} \\
        ${simple_repeat_arg} \\
        ${prefix} \\
        ${tumor_bam} \\
        ${ref}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.nanomonsv.result.txt
    touch ${prefix}.nanomonsv.result.vcf
    touch ${prefix}.nanomonsv.supporting_read.txt
    touch ${prefix}.nanomonsv.sbnd.result.txt
    touch ${prefix}.nanomonsv.sbnd.result.vcf
    """
}
