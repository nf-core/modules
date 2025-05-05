process NEXTCLADE_DATASETGET {
    tag "$dataset"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/936786744b34cf016b948026a6b4e9489011424e15c28dfb2f7d03c31bb4afb5/data' :
        'community.wave.seqera.io/library/nextclade:3.11.0--155203da8341cfe6' }"

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
        nextclade-tag: \$(grep "tag" $dataset/pathogen.json | sed -n 's/.*"tag": "\\([0-9-]\\+Z\\)".*/\\1/p')
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
