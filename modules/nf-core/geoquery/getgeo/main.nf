process GEOQUERY_GETGEO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c367777a9f3947bdfc464c3cec52d33d2862c5e909fa77e13455f8c1d177767b/data' :
        'community.wave.seqera.io/library/bioconductor-geoquery_r-base:8547a00389cd568a' }"

    input:
    tuple val(meta), val(querygse)

    output:
    tuple val(meta), path("*.rds")            , emit: rds
    tuple val(meta), path("*matrix.tsv")      , emit: expression
    tuple val(meta), path("*annotation.tsv")  , emit: annotation
    path "versions.yml"                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'getgeo.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    touch ${prefix}.matrix.tsv
    touch ${prefix}.annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'R.Version()\$version.string' | sed -n 's|\\[1\\] "R version \\(.*\\) (.*|\\1|p')
        bioconductor-geoquery: \$(Rscript -e 'packageVersion("GEOquery")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """
}
