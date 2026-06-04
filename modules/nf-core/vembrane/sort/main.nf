process VEMBRANE_SORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vembrane:2.4.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/vembrane:2.4.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf)
    val expression

    output:
    tuple val(meta), path("*.{vcf,bcf,bcf.gz}"), emit: vcf
    tuple val("${task.process}"), val('vembrane'), eval("vembrane --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_vembrane


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
