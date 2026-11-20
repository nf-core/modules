process CNMF {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/901d99e9e0a403c38effe4db0a0e78fcb8a3516125d6a5c328ac97556c28ea11/data'
        : 'community.wave.seqera.io/library/cnmf:1.7.1--5ebca644e705ff3c'}"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.feature_loadings.csv")  , emit: feature_loadings
    tuple val(meta), path("*.sample_factors.csv")    , emit: sample_factors
    tuple val(meta), path("*.gene_spectra_score.csv"), emit: gene_spectra_score
    tuple val(meta), path("*.stats.csv")             , emit: stats
    tuple val(meta), path("*.merged_spectra.df.npz") , emit: merged_spectra
    tuple val(meta), path("*.clustering.png")        , emit: clustering_plot, optional: true
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cnmf.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feature_loadings.csv
    touch ${prefix}.sample_factors.csv
    touch ${prefix}.gene_spectra_score.csv
    touch ${prefix}.stats.csv
    touch ${prefix}.merged_spectra.df.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnmf: \$(python -c "import importlib.metadata as m; print(m.version('cnmf'))")
        python: \$(python -c "import platform; print(platform.python_version())")
    END_VERSIONS
    """
}
