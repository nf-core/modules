process TAXONKIT_LINEAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.18.0--h9ee0642_0':
        'biocontainers/taxonkit:0.18.0--h9ee0642_0' }"

    input:
    tuple val(meta), val(taxid), path(taxidfile)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert (!taxid && taxidfile) || (taxid && !taxidfile)
    """
    taxonkit \\
        lineage \\
        $args \\
        --data-dir $taxdb \\
        --threads $task.cpus \\
        --out-file ${prefix}.tsv \\
        ${taxid? "<<< '$taxid'": taxidfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """
}
