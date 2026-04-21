process DECOUPLER_DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b5d8e776c3b3a00c25d9a94adbde00008018e5f9eb4aa73c17cb7b886ddc1ae3/data' :
        'community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_python_scanpy:c539733016e15bbc' }"

    input:
    tuple val(meta), path(mat)
    tuple val(meta2), path(net)
    tuple val(meta3), path(annot)

    output:
    tuple val(meta), path("*estimate_decoupler.tsv"), emit: dc_estimate
    tuple val(meta), path("*pvals_decoupler.tsv"), emit: dc_pvals
    tuple val(meta), path("*estimate_decoupler_plot.png"), emit: png
    path "versions.yml", emit: versions_decoupler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'decoupler.py'

    stub:
    """
    touch deseq2_estimate_decoupler.tsv
    touch deseq2_pvals_decoupler.tsv
    touch deseq2_estimate_decoupler_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        decoupler-py: \$(pip show decoupler | sed -n 's/Version: //p')
        matplotlib: \$(pip show matplotlib | sed -n 's/Version: //p')
        pandas: \$(pip show pandas | sed -n 's/Version: //p')
        python: \$(python --version | cut -f 2 -d " ")
        scanpy: \$(pip show scanpy | sed -n 's/Version: //p')
    END_VERSIONS
    """
}
