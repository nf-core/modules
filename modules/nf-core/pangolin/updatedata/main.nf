process PANGOLIN_UPDATEDATA {
    tag "$dbname"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/858e91f6972f0d8d71dae844bf0232656f5d91112b9a5610f559659b33414c86/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_snakemake-minimal:638a1eb68adff9c7' }"

    input:
    val(dbname)

    output:
    path("${prefix}"), emit: db
    tuple val("${task.process}"), val('pangolin-dataset'), eval("pangolin -pv | sed -n 's/pangolin-data \\([0-9.]*\\).*/\\1/p'"), emit: versions_pangolin_dataset, topic: versions
    tuple val("${task.process}"), val('pangolin'), eval("pangolin -v | sed -n 's/pangolin \\([0-9.]*\\).*/\\1/p'"), emit: versions_pangolin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"
    """
    export XDG_CACHE_HOME=/tmp/.cache
    mkdir -p ${prefix}

    pangolin \\
        $args \\
        --update-data \\
        --threads $task.cpus \\
        --datadir ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    mkdir -p ${prefix}
    """
}
