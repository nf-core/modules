process PYDAMAGE_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b5b5289345c54b75f42a54d300d222b9c0de8de5a60ecbe9ad829b6b85f1abd/data'
:         'community.wave.seqera.io/library/pydamage:1.0--1c195c7c48e87ebb' }"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("${prefix}_pydamage_filtered_results.csv"), emit: csv
    tuple val("${task.process}"), val('pydamage'), eval("NUMBA_CACHE_DIR=./tmp pydamage --version | sed -n 's/pydamage, version //p'"), emit: versions_pydamage, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=./tmp
    export MPLCONFIGDIR=./tmp

    pydamage \\
        ${args} \\
        filter \\
        ${csv}

    mv pydamage_results/pydamage_filtered_results.csv ${prefix}_pydamage_filtered_results.csv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pydamage_filtered_results.csv
    """
}
