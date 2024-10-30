process SCVITOOLS_SOLO {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'docker.io/nicotru/scvitools-gpu:cuda-12' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata_scvi-tools:54d2eb2f946e0e43':
        'community.wave.seqera.io/library/anndata_scvi-tools:fa9451a13918eae0' }"

    // Prevent /tmp mount for singularity and set the MATPLOTLIB TMPDIR
    containerOptions  "${workflow.containerEngine == 'singularity' ? '--no-mount tmp --env MPLCONFIGDIR=' + workDir : ''}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl") , emit: predictions
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    batch_key = task.ext.batch_key ?: ""
    max_epochs = task.ext.max_epochs ?: ""
    template 'solo.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvitools: \$(python3 -c 'import scvi; print(scvi.__version__)')
    END_VERSIONS
    """
}
