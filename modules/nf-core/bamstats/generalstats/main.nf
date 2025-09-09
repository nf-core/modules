process BAMSTATS_GENERALSTATS {
    tag "$meta.id"
    label 'process_single'
    conda "bioconda::bamstats=0.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamstats:0.3.5--he881be0_0':
        'biocontainers/bamstats:0.3.5--he881be0_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Several ARGS are available.
    //  -a is a helpfu; one where you can add a BED file
    //  -u, --uniq outputs genomic coverage statistics for uniqely mapped reads
    """
    bamstats \\
        -i $bam \\
        $args \\
        -c $task.cpus \\
        -o ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamstats: \$(echo \$(bamstats --version 2>&1) | sed 's/^.*bamstats == version://; s/Using.*\$//' | sed 's/built.*//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamstats: \$(echo \$(bamstats --version 2>&1) | sed 's/^.*bamstats == version://; s/Using.*\$//' | sed 's/built.*//' )
    END_VERSIONS
    """
}
