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
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "${args}"
    touch ${prefix}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
