process CUSTOM_CLUSTERVISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64297e13d9d4f05ce543656e943a023735e7cb252d7534ccc9134d8b40423083/data' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:826e4ab1361ff931' }"

    input:
    tuple val(meta), path(features), path(clusters)

    output:
    tuple val(meta), path("*.umap.tsv"), emit: umap_tsv
    tuple val(meta), path("*.tsne.tsv"), emit: tsne_tsv
    tuple val(meta), path("*.umap.png"), emit: umap_png, optional: true
    tuple val(meta), path("*.tsne.png"), emit: tsne_png, optional: true
    path "versions.yml"                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cluster_viz.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umap.tsv
    touch ${prefix}.tsne.tsv
    touch ${prefix}.umap.png
    touch ${prefix}.tsne.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c "from importlib.metadata import version; print(version('matplotlib'))")
        numpy: \$(python3 -c "from importlib.metadata import version; print(version('numpy'))")
        pandas: \$(python3 -c "from importlib.metadata import version; print(version('pandas'))")
        scikit-learn: \$(python3 -c "from importlib.metadata import version; print(version('scikit-learn'))")
        seaborn: \$(python3 -c "from importlib.metadata import version; print(version('seaborn'))")
        umap-learn: \$(python3 -c "from importlib.metadata import version; print(version('umap-learn'))")
    END_VERSIONS
    """
}
