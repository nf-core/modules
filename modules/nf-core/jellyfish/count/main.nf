process JELLYFISH_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmer-jellyfish:2.3.1--py310h184ae93_5':
        'quay.io/biocontainers/kmer-jellyfish:2.3.1--py310h184ae93_5' }"

    input:
    tuple val(meta), path(fasta)
    val kmer_length
    val size

    output:
    tuple val(meta), path("${prefix}.jf"), emit: jf
    tuple val("${task.process}"), val("jellyfish"), eval("jellyfish --version |& sed '1!d;s/jellyfish //'"), topic: versions, emit: versions_jellyfish

    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_list = fasta instanceof List ? fasta : [fasta]
    def decompress = fasta_list
        .findAll { input_file -> input_file.getName().endsWith(".gz") }
        .collect { input_file -> "gzip -c -d ${input_file} > ${input_file.getName().replace(".gz", "")}" }
        .join("\n")
    def fasta_names = fasta_list
        .collect { input_file -> input_file.getName().replace(".gz", "") }
        .join(" ")
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${decompress}
    jellyfish \\
        count \\
        $args \\
        -m ${kmer_length} \\
        -s ${size} \\
        -t $task.cpus \\
        -o ${prefix}.jf \\
        ${fasta_names}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.jf
    """
}
