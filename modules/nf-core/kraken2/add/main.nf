process KRAKEN2_ADD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_2' :
        'biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2' }"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names
    path taxonomy_nodes
    path accession2taxid

    output:
    tuple val(meta), path("$prefix"), emit: db
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    cp -L ${taxonomy_names} .
    cp -L ${taxonomy_nodes} .
    cp -L ${accession2taxid} .
    mkdir "${prefix}/taxonomy"
    mv *.{accession2taxid,dmp} "${prefix}/taxonomy"
    kraken2-build \\
        --add-to-library \\
        ${fasta} \\
        --db ${prefix} \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
