process AMPCOMBI2_CLUSTER {
    tag 'ampcombi2'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.2.2--pyhdfd78af_0':
        'biocontainers/ampcombi:0.2.2--pyhdfd78af_0' }"

    input:
    path(summary_file)

    output:
    path("Ampcombi_summary_cluster.tsv")                   , emit: cluster_tsv
    path("Ampcombi_summary_cluster_representative_seq.tsv"), emit: rep_cluster_tsv
    path("Ampcombi_cluster.log")                           , emit: log, optional:true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ampcombi cluster \\
        --ampcombi_summary ${summary_file} \\
        $args \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch Ampcombi_summary_cluster.tsv
    touch Ampcombi_summary_cluster_representative_seq.tsv
    touch Ampcombi_cluster.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
