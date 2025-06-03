process CRABS_IMPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.0.7--pyhdfd78af_0':
        'biocontainers/crabs:1.0.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(accession2taxid)
    tuple val(meta3), path(names)
    tuple val(meta4), path(nodes)
    val(import_format)

    output:
    tuple val(meta), path("*.txt"), emit: crabsdb
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def is_compressed  = fasta.name.endsWith(".gz")
    def fasta_name     = fasta.name.replace(".gz", "")
    def import_fmt_cmd = "--import-format ${import_format}"
    def version_cmd    = "\$(crabs --help 2>/dev/null | grep 'CRABS |' | sed 's/.*CRABS | v\\([0-9.]*\\).*/\\1/')"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        version_cmd    = '1.0.7'
    }
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    crabs --import \\
        --input ${fasta_name} \\
        --output ${prefix}.txt \\
        --acc2tax ${accession2taxid} \\
        --names ${names} \\
        --nodes ${nodes} \\
        ${import_fmt_cmd} \\
        ${args}

    rm ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: ${version_cmd}
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def version_cmd  = "\$(crabs --help 2>/dev/null | grep 'CRABS |' | sed 's/.*CRABS | v\\([0-9.]*\\).*/\\1/')"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        version_cmd  = '1.0.7'
    }
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: ${version_cmd}
    END_VERSIONS
    """
}
