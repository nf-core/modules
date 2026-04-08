process DIAMOND_MAKEDB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.16--h13889ed_0'
        : 'biocontainers/diamond:2.1.16--h13889ed_0'}"

    input:
    tuple val(meta), path(fasta)
    path taxonmap
    path taxonnodes
    path taxonnames

    output:
    tuple val(meta), path("*.dmnd"), emit: db
    tuple val("${task.process}"), val('diamond'), eval("diamond --version | sed 's/diamond version //g'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def insert_taxonmap = taxonmap ? "--taxonmap ${taxonmap}" : ""
    def insert_taxonnodes = taxonnodes ? "--taxonnodes ${taxonnodes}" : ""
    def insert_taxonnames = taxonnames ? "--taxonnames ${taxonnames}" : ""

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    diamond \\
        makedb \\
        --threads ${task.cpus} \\
        --in  ${fasta_name} \\
        -d ${prefix} \\
        ${args} \\
        ${insert_taxonmap} \\
        ${insert_taxonnodes} \\
        ${insert_taxonnames}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "${args}"
    touch ${prefix}.dmnd
    """
}
