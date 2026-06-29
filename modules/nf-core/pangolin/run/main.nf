process PANGOLIN_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f10f9bd63d24e2382406b2b348c65ae1ac74118e6eb17f5e30b310fbc1bc4b9/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_pip_snakemake:7987089a2a044e70' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.csv'), emit: report
    tuple val("${task.process}"), val('pangolin'), eval("pangolin -v | sed -n 's/pangolin \\([0-9.]*\\).*/\\1/p'"), emit: versions_pangolin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_command = db ? "--datadir ${db}" : ''
    """
    export XDG_CACHE_HOME=/tmp/.cache

    pangolin \\
        $fasta\\
        $db_command \\
        --outfile ${prefix}.pangolin.csv \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=/tmp/.cache

    touch ${prefix}.pangolin.csv
    """
}
