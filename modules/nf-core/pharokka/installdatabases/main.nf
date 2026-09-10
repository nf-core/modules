process PHAROKKA_INSTALLDATABASES {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.9.1--pyhdfd78af_0':
        'quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_0' }"

    output:
    path("${prefix}/"), emit: pharokka_db
    tuple val("${task.process}"), val('pharokka'), eval("pharokka.py --version"), topic: versions, emit: versions_pharokka

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'pharokka_db'
    """
    install_databases.py \\
        --outdir ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'pharokka_db'
    """
    mkdir -p ${prefix}
    touch ${prefix}/1Aug2023_data.tsv
    touch ${prefix}/1Aug2023_genomes.fa.msh
    touch ${prefix}/CARD{,.dbtype,.index,.lookup,.source}
    touch ${prefix}/CARD_h{,.dbtype,.index}
    echo "" | gzip > ${prefix}/VFDB_setB_pro.fas.gz
    touch ${prefix}/VFDBclusterRes_cluster.tsv
    touch ${prefix}/VFDBclusterRes_rep_seq.fasta
    touch ${prefix}/all_phrogs.h3m
    touch ${prefix}/aro_index.tsv
    touch ${prefix}/phrog_annot_v4.tsv
    touch ${prefix}/phrogs_db{,.dbtype,.index}
    touch ${prefix}/phrogs_profile_db{,.dbtype,.index}
    touch ${prefix}/phrogs_profile_db_consensus{,.dbtype,.index}
    touch ${prefix}/phrogs_profile_db_h{,.index}
    touch ${prefix}/phrogs_profile_db_seq{,.dbtype,.index}
    touch ${prefix}/phrogs_profile_db_seq_h{,.index}
    touch ${prefix}/vfdb{,.dbtype,.index,.lookup,.source}
    touch ${prefix}/vfdb_h{,.dbtype,.index}
    """
}
