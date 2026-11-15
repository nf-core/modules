process SVDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdss:2.1.1--he17396a_0' :
        'quay.io/biocontainers/svdss:2.1.1--he17396a_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(sfs)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    SVDSS call \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        --bam ${bam} \\
        --sfs ${sfs} \\
        ${args} \\
        | bgzip -@ ${task.cpus} -c > ${prefix}.vcf.gz
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    echo "" | gzip > ${prefix}.vcf.gz
    """
}
