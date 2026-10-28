process COLABFOLD_DOWNLOADMODELS {
    tag "${model_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2de6276b0c23cc1aaab4af2aa1336fa943b0c45d43ce4ed483911fb99a3db33/data'
        : 'community.wave.seqera.io/library/colabfold_dm-haiku_flax:50a22b7e73876edb'}"

    input:
    val(model_type)

    output:
    tuple val(meta), path("weights"), emit: weights
    tuple val("${task.process}"), val('colabfold'), eval("python -c 'import importlib.metadata as m; print(m.version(\"colabfold\"))'"), topic: versions, emit: versions_colabfold

    when:
    task.ext.when == null || task.ext.when

    script:
    meta = [ id: model_type ]
    """
    # Same downloader colabfold_batch calls itself; it also writes the
    # params/download*_finished.txt marker that colabfold_batch checks for
    python -c "from pathlib import Path; from colabfold.download import download_alphafold_params; download_alphafold_params('${model_type}', Path('weights'))"
    """

    stub:
    meta = [ id: model_type ]
    def markers = [
        alphafold2_multimer_v1: 'download_complexes_multimer_v1_finished.txt',
        alphafold2_multimer_v2: 'download_complexes_multimer_v2_finished.txt',
        alphafold2_multimer_v3: 'download_complexes_multimer_v3_finished.txt',
        deepfold_v1           : 'download_deepfold-v1_finished.txt',
    ]
    """
    mkdir -p weights/params
    touch weights/params/${markers[model_type] ?: 'download_finished.txt'}
    """
}
