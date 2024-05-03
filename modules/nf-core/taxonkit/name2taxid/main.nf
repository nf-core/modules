process TAXONKIT_NAME2TAXID {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.15.1--h9ee0642_0':
        'biocontainers/taxonkit:0.15.1--h9ee0642_0' }"

    input:
    tuple val(meta), val(name), path(names_txt)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert (!name && names_txt) || (name && !names_txt)
    """
    taxonkit \\
        name2taxid \\
        $args \\
        --data-dir $taxdb \\
        --threads $task.cpus \\
        --out-file ${prefix}.tsv \\
        ${name? "<<< '$name'": names_txt}

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
