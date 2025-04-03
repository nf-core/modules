process TRYCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trycycler:0.5.3--pyhdfd78af_0':
        'biocontainers/trycycler:0.5.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*/*.fastq.gz") , emit: subreads
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    trycycler \\
        subsample \\
        $args \\
        --reads ${reads} \\
        --threads $task.cpus \\
        --out_dir ${prefix}

    gzip $args2 ${prefix}/*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(trycycler --version | sed 's/Trycycler v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Pull out the count integer from the args string, if provided
    // Otherwise use default value of 12
    def matches = (args =~ /--count\s+(\d+)/)
    def count = matches ? matches[0][1] as Integer : 12
    """
    for n in \$(seq $count); do
        printf -v i "%02d" \$n
        mkdir -p ${prefix}
        echo "" | gzip > ${prefix}/sample_\${i}.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(trycycler --version | sed 's/Trycycler v//')
    END_VERSIONS
    """
}
