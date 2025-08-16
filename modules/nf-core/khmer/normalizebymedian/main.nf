process KHMER_NORMALIZEBYMEDIAN {
    tag "${meta.id}"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(fastq_paired), path(fastq_unpaired)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_normalized"
    assert "${prefix}.fastq.gz" != fastq_paired.name   : "KHMER_NORMALIZEBYMEDIAN: The output file name must be different from the input file name. Please set a different prefix."
    assert "${prefix}.fastq.gz" != fastq_unpaired.name : "KHMER_NORMALIZEBYMEDIAN: The output file name must be different from the input file name. Please set a different prefix."
    def paired   = fastq_paired   ? "--paired"                     : ""
    def unpaired = fastq_unpaired ? "--unpaired ${fastq_unpaired}" : ""
    def input_sequence_filename = fastq_paired ? fastq_paired : fastq_unpaired
    """
    normalize-by-median.py \\
        -M ${task.memory.toGiga()}e9 \\
        --gzip $args \\
        -o ${prefix}.fastq.gz \\
        ${paired} \\
        ${unpaired} \\
        ${input_sequence_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_normalized"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
