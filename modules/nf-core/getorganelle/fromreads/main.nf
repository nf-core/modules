process GETORGANELLE_FROMREADS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::getorganelle=1.7.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(organelle_type), path(db)  // getOrganelle has a database and config file

    output:
    tuple val(meta), path("results/${prefix}.${organelle_type}.fasta.gz"), emit: fasta, optional: true
    path "results/*"                                                     , emit: etc // the rest of the result files
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    get_organelle_from_reads.py \\
        $args \\
        --prefix ${meta.id}. \\
        -F $organelle_type \\
        --config-dir $db \\
        -t $task.cpus \\
        -1 ${fastq[0]} \\
        -2 ${fastq[1]} \\
        -o results

    if [ -f "\$(find results -name ${prefix}.${organelle_type}*graph1.1*fasta )" ]; then
        cp results/${prefix}.${organelle_type}*graph1.1*fasta results/${prefix}.${organelle_type}.fasta
        gzip results/${prefix}.${organelle_type}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_from_reads.py --version | sed 's/^GetOrganelle v//g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > results/${prefix}.${organelle_type}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_from_reads.py --version | sed 's/^GetOrganelle v//g' )
    END_VERSIONS
    """
}
