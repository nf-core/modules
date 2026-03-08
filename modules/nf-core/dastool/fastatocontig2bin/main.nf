process DASTOOL_FASTATOCONTIG2BIN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/das_tool:1.1.7--r44hdfd78af_1'
        : 'biocontainers/das_tool:1.1.7--r44hdfd78af_1'}"

    input:
    tuple val(meta), path(fasta)
    val extension

    output:
    tuple val(meta), path("*.tsv"), emit: fastatocontig2bin
    tuple val("${task.process}"), val('dastool'), eval('DAS_Tool --version 2>&1 | grep "DAS Tool" | sed "s/DAS Tool //"'), topic: versions, emit: versions_dastool

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = extension ? extension : "fasta"
    def clean_fasta = fasta.toString() - ".gz"
    def decompress_fasta = fasta.toString() == clean_fasta ? "" : "gunzip -q -f ${fasta}"
    """
    ${decompress_fasta}

    Fasta_to_Contig2Bin.sh \\
        ${args} \\
        -i . \\
        -e ${file_extension} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
