process BARRNAP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::barrnap=0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/barrnap:0.9--hdfd78af_4':
        'biocontainers/barrnap:0.9--hdfd78af_4' }"

    input:
    tuple val(meta), path(reads), val(dbname)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    db         = dbname ? "${dbname}" : 'bac'

    """
    barrnap \\
        $args \\
         --threads $task.cpus \\
        --kingdom $db \\
        $reads 
        > rrna_${db}.gff


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(barrnap --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(barrnap --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
