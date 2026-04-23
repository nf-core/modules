process VCF2CIRCOS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vcf2circos:1.2.0--pyhdfd78af_0'
        : 'biocontainers/vcf2circos:1.2.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(vcf_index), val(output_extension)
    tuple val(meta2), path(vcf2circos_config)

    output:
    tuple val(meta), path("*.{html,png,jpg,jpeg,webp,svg,pdf,eps,json}"), emit: circos
    tuple val("${task.process}"), val('vcf2circos'), eval("vcf2circos --help | sed -n 's/^Version: //p'"), topic: versions, emit: versions_vcf2circos

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2circos \\
        ${args} \\
        --options ${vcf2circos_config} \\
        --input ${vcf} \\
        -o ${prefix}.${output_extension}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.${output_extension}
    """
}
