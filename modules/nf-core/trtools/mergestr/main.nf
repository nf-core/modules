process TRTOOLS_MERGESTR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trtools:6.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/trtools:6.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('trtools'), eval("mergeSTR --version | sed 's/mergeSTR //'"), topic: versions, emit: versions_trtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_mergestr"
    def input = vcfs.sort { vcf -> vcf.toString() }.join(",")

    """
    mergeSTR \\
        --vcfs ${input} \\
        --out ${prefix} \\
        ${args}

    bgzip -f ${prefix}.vcf
    tabix -f -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
