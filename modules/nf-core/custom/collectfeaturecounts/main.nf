process CUSTOM_COLLECTFEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/143121bf17b3f4e9539d8cd17aaaa428aae54ce827471c4fc6399465263efa2e/data' :
        'community.wave.seqera.io/library/r-base_r-dplyr_r-readr_r-purrr_pruned:0f879b99d6a89834' }"

    input:
    tuple val(meta), path(inputfiles, stageAs: "input/*")

    output:
    tuple val(meta), path("${prefix}.counts.tsv.gz"), emit: counts
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'collectfeaturecounts.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    gzip ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
        stringr: 1.4.0
        dtplyr: 1.1.0
        data.table: 1.14.0
    END_VERSIONS
    """
}
