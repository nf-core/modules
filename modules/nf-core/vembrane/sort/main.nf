process VEMBRANE_SORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a1bc5a7c85eebca118ae436bda77e005cc6c1ca0e60687185fb1d2f195e5eaf/data'
        : 'community.wave.seqera.io/library/vembrane:2.6.2--0bc429af0e13a402'}"

    input:
    tuple val(meta), path(vcf)
    val expression

    output:
    tuple val(meta), path("*.{vcf,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('vembrane'), eval("vembrane --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_vembrane

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains('--output-fmt bcf')
        ? 'bcf'
        : args.contains('--output-fmt bcf.gz') ? 'bcf.gz' : 'vcf'
    """
    vembrane sort \\
        ${args} \\
        --output ${prefix}.${suffix} \\
        '${expression}' \\
        ${vcf}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains('--output-fmt bcf')
        ? 'bcf'
        : args.contains('--output-fmt bcf.gz') ? 'bcf.gz' : 'vcf'
    """
    touch ${prefix}.${suffix}
    """
}
