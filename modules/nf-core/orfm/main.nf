process ORFM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ef/efb1c20e288683f9dbbe6240e81869afb5339b41a63706e5fe573658c5c7621f/data':
        'community.wave.seqera.io/library/orfm:2.1.1--d2f1f66d5b27e4e9' }"

    input:
    tuple val(meta), path(sequence)
    tuple val(meta2), val(min_len), val(codon_table_id)


    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('orfm'), eval("orfm --version | sed 's/orfm //'"), topic: versions, emit: versions_orfm

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def codon_table = codon_table_id ?: '1'
    // check if the provided codon table is a valid integer
    if (codon_table.toInteger() <= 0 || codon_table.toInteger() > 26) {
        error "The provided codon table (${codon_table}) should be a valid integer between 1 and 26."
    }

    def min_orflen = min_len ?: '96'
    // check if the provided minimum ORF length is a multiple of 3
    if (min_orflen.toInteger() % 3 != 0) {
        error "The provided minimum ORF length (${min_orflen}) should be in a multiple of 3."
    }

    """
    orfm \\
        ${sequence} \\
        -m ${min_orflen} \\
        -c ${codon_table} \\
        -t ${prefix}.fasta \\
        -s \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """    
    touch ${prefix}.fasta
    """
}
