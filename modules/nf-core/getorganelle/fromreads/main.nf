process GETORGANELLE_FROMREADS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fastq)
    tuple val(organelle_type), path(db)

    output:
    tuple val(meta), path("results/${prefix}.${organelle_type}.fasta.gz"), emit: fasta, optional: true
    tuple val(meta), path("results/*"), emit: etc
    tuple val("${task.process}"), val('getorganelle'), eval("get_organelle_from_reads.py --version |& sed 's/^GetOrganelle v//'"), topic: versions, emit: versions_getorganelle

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    get_organelle_from_reads.py \\
        ${args} \\
        --prefix ${meta.id}. \\
        -F ${organelle_type} \\
        --config-dir ${db} \\
        -t ${task.cpus} \\
        -1 ${fastq[0]} \\
        -2 ${fastq[1]} \\
        -o results

    if [ -f "\$(find results -name ${prefix}.${organelle_type}*graph1.1*fasta )" ]; then
        cp results/${prefix}.${organelle_type}*graph1.1*fasta results/${prefix}.${organelle_type}.fasta
        gzip results/${prefix}.${organelle_type}.fasta
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '${args}'
    mkdir results
    touch results/test_1.fastq
    touch results/test_2.fastq
    echo "" | gzip > results/${prefix}.${organelle_type}.fasta.gz
    """
}
