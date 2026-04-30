process HAPLOCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplocheck:1.3.3--h4a94de4_0':
        'quay.io/biocontainers/haplocheck:1.3.3--h4a94de4_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.txt") , emit: txt
    tuple val(meta), path("*.html"), emit: html
    tuple val("${task.process}"), val('haplocheck'), eval("haplocheck --version 2>&1 | head -1 | awk '{print \$2}'"), emit: versions_haplocheck, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    haplocheck --raw --out $prefix $vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.txt
    touch ${prefix}.html
    """
}
