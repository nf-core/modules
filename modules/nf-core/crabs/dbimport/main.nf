process CRABS_DBIMPORT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:0.1.1--pyhb7b1952_0':
        'biocontainers/crabs:0.1.1--pyhb7b1952_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name    = fasta.name.replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    crabs db_import \\
        --input ${fasta_name} \\
        --output ${prefix}.crabsdb.fa \\
        $args

    rm ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '0.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: $VERSION
    END_VERSIONS
    """
}
