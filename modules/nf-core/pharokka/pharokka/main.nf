process PHAROKKA_PHAROKKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.5.1--pyhdfd78af_0 ':
        'biocontainers/pharokka:1.5.1--pyhdfd78af_0 ' }"

    input:
    tuple val(meta), path(phage_fasta)
    path pharokka_db

    output:
    tuple val(meta), path("*/*.gbk")                                , emit: genbank
    tuple val(meta), path("*/*.log")                                , emit: log
    tuple val(meta), path("*/*_cds_functions.tsv")                  , emit: cds_functions
    tuple val(meta), path("*/*top_hits_card.tsv")                   , emit: card
    tuple val(meta), path("*/*top_hits_vfdb.tsv")                   , emit: vfdb
    tuple val(meta), path("*/*_top_hits_mash_inphared.tsv")         , emit: mash            , optional: true
    tuple val(meta), path("*/*_genome_terminase_reoriented.fasta")  , emit: reoriented      , optional: true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pharokka.py \\
        --infile ${phage_fasta} \\
        --outdir ${meta.id}_pharokka \\
        --database ${pharokka_db} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${meta.id}_pharokka
    touch ${meta.id}_pharokka/${meta.id}.gbk
    touch ${meta.id}_pharokka/${meta.id}.log
    touch ${meta.id}_pharokka/${meta.id}_cds_functions.tsv
    touch ${meta.id}_pharokka/${meta.id}_top_hits_card.tsv
    touch ${meta.id}_pharokka/top_hits_vfdb.tsv
    touch ${meta.id}_pharokka/${meta.id}_top_hits_inphared
    touch ${meta.id}_pharokka/${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """
}
