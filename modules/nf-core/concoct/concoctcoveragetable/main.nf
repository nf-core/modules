process CONCOCT_CONCOCTCOVERAGETABLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py39h8907335_8'
        : 'biocontainers/concoct:1.1.0--py39h8907335_8'}"

    input:
    tuple val(meta), path(bed), path(bamfiles), path(baifiles)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('concoct'), eval('concoct --version | cut -d " " -f2'), emit: versions_concoct, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    concoct_coverage_table.py \\
        ${args} \\
        ${bed} \\
        ${bamfiles} \\
        > ${prefix}.tsv

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    """
}
