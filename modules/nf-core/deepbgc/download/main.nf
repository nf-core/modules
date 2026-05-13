process DEEPBGC_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    output:
    path "deepbgc_db/", emit: db
    tuple val("${task.process}"), val('deepbgc'), eval("deepbgc info 2>&1 | sed '6!d;s/.*= version //;s/ .*//'"), emit: versions_deepbgc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    export DEEPBGC_DOWNLOADS_DIR='./deepbgc_db'

    deepbgc \\
        download
    """

    stub:
    """
    mkdir -p deepbgc_db/0.1.0/classifier deepbgc_db/0.1.0/detector deepbgc_db/common

    touch deepbgc_db/common/Pfam-A.31.0.clans.tsv
    touch deepbgc_db/common/Pfam-A.31.0.hmm.{,h3f,h3i,h3m,h3p}

    touch deepbgc_db/0.1.0/classifier/product_{activity,class}.pkl
    touch deepbgc_db/0.1.0/detector/clusterfinder_{geneborder,original,retrained}.pkl
    touch deepbgc_db/0.1.0/detector/deepbgc.pkl
    """
}
