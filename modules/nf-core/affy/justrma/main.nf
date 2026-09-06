process AFFY_JUSTRMA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/62/62d22bc460807a1a4ded40e5b7a391aa6f2dac189d4153d684472d65333ca8d4/data':
        'community.wave.seqera.io/library/bioconductor-affy_r-base:dd8a5ecd6fc301b3' }"

    input:
    tuple val(meta), path(samplesheet), path(celfiles_dir)
    tuple val(meta2), path(description)

    output:
    tuple val(meta), path("*.rds")            , emit: rds
    tuple val(meta), path("*matrix.tsv")      , emit: expression
    tuple val(meta), path("*.annotation.tsv") , emit: annotation, optional: true
    path "versions.yml", emit: versions_affy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'affy_justrma.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eset.rds
    touch ${prefix}_matrix.tsv
    touch R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        bioconductor-affy: \$(Rscript -e "cat(as.character(packageVersion('affy')))")
    END_VERSIONS
    """
}
