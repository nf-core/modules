process DREP_COMPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/drep:3.6.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/drep:3.6.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${prefix}"), emit: directory
    tuple val("${task.process}"), val("drep"), eval("dRep | sed '2!d;s/.*v//g;s/ .*//g'"), emit:versions_drep, topic:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    dRep \\
        compare \\
        ${prefix} \\
        -p ${task.cpus} \\
        ${args} \\
        -g fastas/*
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p ${prefix}/data ${prefix}/data_tables ${prefix}/dereplicated_genomes ${prefix}/figures ${prefix}/log
    mkdir -p ${prefix}/data/Clustering_files ${prefix}/data/fastANI_files ${prefix}/data/MASH_files

    touch ${prefix}/data/Clustering_files/primary_linkage.pickle
    touch ${prefix}/data/fastANI_files/fastANI_out_{jfkgcosewp,szxnawefbt}{,.matrix}
    touch ${prefix}/data/MASH_files/chunk_0_MASH_table.tsv

    touch ${prefix}/data_tables/{Bdb,Cdb,Mdb,Ndb}.csv

    touch ${prefix}/figures/{Clustering_scatterplots,Primary_clustering_dendrogram}.pdf

    touch ${prefix}/log/cluster_arguments.json
    touch ${prefix}/log/logger.log
    touch ${prefix}/log/warnings.txt
    """
}
