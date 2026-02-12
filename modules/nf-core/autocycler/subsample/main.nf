
process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

    input:
    tuple val(meta), path(reads)
    val genome_size

    output:
    tuple val(meta), path("$prefix/*.fastq.gz"), emit: subsampled_reads
    tuple val("${task.process}"), val("autocycler"), eval("autocycler --version |  sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // fix random seed for reproducibility if not specified in command line
    if (!(args ==~ /.*--seed.*/)) {args += " --seed 42"}
    """
    autocycler subsample \\
        $args \\
        --reads $reads \\
        --out_dir ${prefix} \\
        --genome_size $genome_size

    gzip $prefix/*.fastq
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """

    mkdir $prefix
    echo | gzip > ${prefix}/sample_00.fastq.gz
    """
}
