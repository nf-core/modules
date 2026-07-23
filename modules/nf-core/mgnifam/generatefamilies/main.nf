process MGNIFAM_GENERATEFAMILIES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/aa/aa77afddce5309d57c80ddc131bb6ca7232d1227abe608f87011cb1caf36c4ff/data' :
        'community.wave.seqera.io/library/pip_mgnifam:b23278d764db6a8f' }"

    input:
    tuple val(meta), path(clusters_chunk), path(fasta_file)

    output:
    tuple val(meta), path("${prefix}/seed_msa/*.sto.gz")       , emit: seed_msa  , optional: true
    tuple val(meta), path("${prefix}/full_msa/*.sto.gz")       , emit: full_msa  , optional: true
    tuple val(meta), path("${prefix}/hmm/*.hmm.gz")            , emit: hmm       , optional: true
    tuple val(meta), path("${prefix}/rf/*.txt")                , emit: rf        , optional: true
    tuple val(meta), path("${prefix}/${prefix}_families.tsv")  , emit: tsv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_metadata.csv")  , emit: csv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}.log")           , emit: log       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_reps.fasta.gz") , emit: reps_fasta, optional: true
    tuple val(meta), path("${prefix}/${prefix}_successful.txt"), emit: successful, optional: true
    tuple val(meta), path("${prefix}/${prefix}_discarded.csv") , emit: discarded , optional: true
    tuple val(meta), path("${prefix}/${prefix}_converged.txt") , emit: converged , optional: true
    tuple val("${task.process}"), val('mgnifam'), eval("mgnifam --version 2>&1"), topic: versions, emit: versions_mgnifam
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mgnifam generate_families \\
        --clusters_chunk ${clusters_chunk} \\
        --fasta_file ${fasta_file} \\
        --output_dir ${prefix} \\
        --cpus ${task.cpus} \\
        --chunk_num ${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}/seed_msa ${prefix}/full_msa ${prefix}/hmm ${prefix}/rf
    python3 -c "import gzip; gzip.open('${prefix}/seed_msa/${prefix}_1.sto.gz', 'wb').close()"
    python3 -c "import gzip; gzip.open('${prefix}/full_msa/${prefix}_1.sto.gz', 'wb').close()"
    python3 -c "import gzip; gzip.open('${prefix}/hmm/${prefix}_1.hmm.gz', 'wb').close()"
    touch ${prefix}/rf/${prefix}_1.txt
    touch ${prefix}/${prefix}_families.tsv
    touch ${prefix}/${prefix}_metadata.csv
    touch ${prefix}/${prefix}.log
    python3 -c "import gzip; gzip.open('${prefix}/${prefix}_reps.fasta.gz', 'wb').close()"
    touch ${prefix}/${prefix}_successful.txt
    touch ${prefix}/${prefix}_discarded.csv
    touch ${prefix}/${prefix}_converged.txt
    """
}
