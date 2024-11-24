process HOMER_FINDPEAKS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/898620b4f7461c66fc41f8532d3ae1cb6ae9f44d304a59c988893244dbcbe389/data' :
        'community.wave.seqera.io/library/homer_wget_samtools_r-essentials_pruned:60f842302c4fb688' }"

    input:
    tuple val(meta), path(tagDir)
    path uniqmap // optional

    output:
    tuple val(meta), path("*.peaks.txt"), emit: txt
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: uniqmap ? "${meta.id}-${uniqmap.baseName}" : "${meta.id}"
    def uniqmap_flag = uniqmap ? "-uniqmap $uniqmap" : ""
    def VERSION = '4.11' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    findPeaks \\
        $tagDir \\
        $args \\
        -o ${prefix}.peaks.txt \\
        $uniqmap_flag

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.peaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
