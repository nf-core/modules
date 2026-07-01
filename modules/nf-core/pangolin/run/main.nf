process PANGOLIN_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/858e91f6972f0d8d71dae844bf0232656f5d91112b9a5610f559659b33414c86/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_snakemake-minimal:638a1eb68adff9c7' }"

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
