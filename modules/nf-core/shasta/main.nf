process SHASTA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::shasta=0.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shasta:0.8.0--h7d875b9_0':
        'biocontainers/shasta:0.8.0--h7d875b9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_Assembly.fasta.gz"), emit: assembly
    tuple val(meta), path("*_Assembly.gfa.gz")  , emit: gfa
    tuple val(meta), path("ShastaRun/")                 , emit: results
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model  = "${meta.model}" ?: 'Nanopore-Oct2021'
    """
    # shasta requires uncompressed
    zcat -f $reads > reads.fq

    # run shasta
    shasta \\
        --input reads.fq \\
        --config $model \\
        $args \\
        --threads $task.cpus

    # compress results
    gzip -c ShastaRun/Assembly.fasta > ${prefix}_Assembly.fasta.gz
    gzip -c ShastaRun/Assembly.gfa   > ${prefix}_Assembly.gfa.gz

    # cleanup temp files
    rm reads.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shasta: \$(shasta --version | head -n 1 | cut -f 3 -d " ")
    END_VERSIONS
    """
}
