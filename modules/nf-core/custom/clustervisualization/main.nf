process CUSTOM_CLUSTERVISUALIZATION {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'docker://community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn_pruned:054f91aaa56bd7d5' :
    'community.wave.seqera.io/library/matplotlib_pandas_python_scikit-learn_pruned:054f91aaa56bd7d5' }"
    input:
    tuple val(meta), path(features), path(clusters)

    output:
    tuple val(meta), path("*.umap.tsv") , emit: umap_tsv
    tuple val(meta), path("*.tsne.tsv") , emit: tsne_tsv
    tuple val(meta), path("*.umap.png") , emit: umap_png, optional: true
    tuple val(meta), path("*.tsne.png") , emit: tsne_png, optional: true
    path "versions.yml"                 , emit: versions, topic: versions

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
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python3 -c "import seaborn; print(seaborn.__version__)")
        umap-learn: \$(python3 -c "import umap; print(umap.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
