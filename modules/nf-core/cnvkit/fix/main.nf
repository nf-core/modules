process CNVKIT_FIX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.13--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(targetcoverage), path(antitargetcoverage), path(reference)

    output:
    tuple val(meta), path("*.cnr"), emit: copy_number_ratio

    //version
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cnvkit.py fix \\
        ${targetcoverage} \\
        ${antitargetcoverage} \\
        ${reference} \\
        ${args} \\
        --output ${prefix}.cnr
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.cnr
    """
}
