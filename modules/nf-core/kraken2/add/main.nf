process KRAKEN2_ADD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f8c4015c836dd3ce5c118cfed97ec8259bab9e9d:2e0b144854b4a3d69b5df7a0340a60db846cc8bf-0' :
        'biocontainers/mulled-v2-f8c4015c836dd3ce5c118cfed97ec8259bab9e9d:2e0b144854b4a3d69b5df7a0340a60db846cc8bf-0' }"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names, stageAs: 'taxonomy/*'
    path taxonomy_nodes, stageAs: 'taxonomy/*'
    path accession2taxid, stageAs: 'taxonomy/*'

    output:
    tuple val(meta), path("$prefix"), emit: db
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    mv "taxonomy" ${prefix}

    echo ${fasta} |\\
    tr -s " " "\\012" |\\
    xargs -I {} -n1 kraken2-build \\
        --add-to-library {} \\
        --db ${prefix} \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

}
