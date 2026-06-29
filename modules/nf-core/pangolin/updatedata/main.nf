process PANGOLIN_UPDATEDATA {
    tag "$dbname"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f10f9bd63d24e2382406b2b348c65ae1ac74118e6eb17f5e30b310fbc1bc4b9/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_pip_snakemake:7987089a2a044e70' }"

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
