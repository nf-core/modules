process BARRNAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/barrnap:0.9--hdfd78af_4':
        'biocontainers/barrnap:0.9--hdfd78af_4' }"

    input:
    tuple val(meta), path(fasta), val(dbname)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val("${task.process}"), val('barrnap'), eval('barrnap --version 2>&1 | sed "s/barrnap //g"'), emit: versions_barrnap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    db         = dbname ? "${dbname}" : 'bac'
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    $gunzip

    barrnap \\
        $args \\
        --threads $task.cpus \\
        --kingdom $db \\
        $input \\
        > ${prefix}_${db}.gff

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    db = dbname ? "${dbname}" : 'bac'
    """
    touch ${prefix}_${db}.gff
    """
}
