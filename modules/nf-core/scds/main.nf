process SCDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f34b3007ea36d52c6cb43692226ab05e79e837b8400693e665214c8f4e2708c/data' :
        'community.wave.seqera.io/library/bioconductor-scds_xgboost:ededce7b92e37374' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: predictions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${rds}" == "${prefix}.rds") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    template 'scds.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${rds}" == "${prefix}.rds") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.rds
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | grep -oP "\\d+\\.\\d+\\.\\d+")
        scds: \$(Rscript -e "library(scds); cat(as.character(packageVersion('scds')), '\\n', sep='')")
        singlecellexperiment: \$(Rscript -e "library(SingleCellExperiment); cat(as.character(packageVersion('SingleCellExperiment')), '\\n', sep='')")
    END_VERSIONS
    """
}
