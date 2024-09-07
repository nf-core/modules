process NEXTCLADE_DATASETGET {
    tag "$dataset"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:3.8.2--h9ee0642_0' :
        'biocontainers/nextclade:3.8.2--h9ee0642_0' }"

    input:
    val dataset
    val tag

    output:
    path "$prefix"     , emit: dataset
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}"
    def version = tag ? "--tag ${tag}" : ''
    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        $version \\
        --output-dir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${dataset}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/CHANGELOG.md
    touch ${prefix}/README.md
    touch ${prefix}/genome_annotation.gff3
    touch ${prefix}/pathogen.json
    touch ${prefix}/reference.fasta
    touch ${prefix}/sequences.fasta
    touch ${prefix}/tree.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """

}
