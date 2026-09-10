process CUSTOM_COLLECTFEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bfe4a872ba15a90cc12fb24aa96ddda852d5b295c214684d9c8fdf3dc02af535/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:e289d008f8e006c5' }"

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
        R: \$(Rscript -e 'cat(paste0(R.Version()[c("major","minor")], collapse = "."))')
        dplyr: \$(Rscript -e 'cat(as.character(packageVersion("dplyr")))')
        readr: \$(Rscript -e 'cat(as.character(packageVersion("readr")))')
        stringr: \$(Rscript -e 'cat(as.character(packageVersion("stringr")))')
        dtplyr: \$(Rscript -e 'cat(as.character(packageVersion("dtplyr")))')
        data.table: \$(Rscript -e 'cat(as.character(packageVersion("data.table")))')
    END_VERSIONS
    """
}
