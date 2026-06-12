process OCTOPUSV_SVCF2VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/octopusv:0.3.2--pyhdfd78af_0' :
        'quay.io/biocontainers/octopusv:0.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(svcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    octopusv svcf2vcf \\
        -i ${svcf} \\
        -o ${prefix}.vcf \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.vcf
    """
}
