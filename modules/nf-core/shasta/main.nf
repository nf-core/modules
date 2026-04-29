process SHASTA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shasta:0.8.0--h7d875b9_0':
        'quay.io/biocontainers/shasta:0.8.0--h7d875b9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_Assembly.fasta.gz"), emit: assembly
    tuple val(meta), path("*_Assembly.gfa.gz")  , emit: gfa
    tuple val(meta), path("ShastaRun/")                 , emit: results
    tuple val("${task.process}"), val('shasta'), eval('shasta --version | head -n 1 | cut -f 3 -d " "'), emit: versions_shasta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: '--config Nanopore-Oct2021'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # shasta requires uncompressed
    zcat -f $reads > reads.fq

    # run shasta
    shasta \\
        --input reads.fq \\
        $args2 \\
        $args \\
        --threads $task.cpus

    # compress results
    gzip -c ShastaRun/Assembly.fasta > ${prefix}_Assembly.fasta.gz
    gzip -c ShastaRun/Assembly.gfa   > ${prefix}_Assembly.gfa.gz

    # cleanup temp files
    rm reads.fq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_Assembly.fasta.gz
    echo "" | gzip > ${prefix}_Assembly.gfa.gz
    mkdir -p ShastaRun
    """
}
