process PANGOLIN_UPDATEDATA {
    tag "$dbname"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f10f9bd63d24e2382406b2b348c65ae1ac74118e6eb17f5e30b310fbc1bc4b9/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_pip_snakemake:7987089a2a044e70' }"

    input:
    val(dbname)

    output:
    path("${prefix}")   , emit: db
    path "versions.yml" , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """
}
