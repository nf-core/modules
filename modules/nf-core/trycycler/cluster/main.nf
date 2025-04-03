process TRYCYCLER_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trycycler:0.5.3--pyhdfd78af_0':
        'biocontainers/trycycler:0.5.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(reads)

    output:
    tuple val(meta), path("*") , emit: cluster_dir
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    trycycler \\
        cluster \\
        $args \\
        --assemblies ${contigs} \\
        --reads ${reads} \\
        --threads $task.cpus \\
        --out_dir ${prefix}

    gzip $args2 ${prefix}/cluster_*/*/*.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(trycycler --version | sed 's/Trycycler v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/cluster_001/1_contigs
    echo "" | gzip > ${prefix}/cluster_001/1_contigs/A_contig_2a.fasta.gz
    touch ${prefix}/contigs.newick
    touch ${prefix}/contigs.phylip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(trycycler --version | sed 's/Trycycler v//')
    END_VERSIONS
    """
}
