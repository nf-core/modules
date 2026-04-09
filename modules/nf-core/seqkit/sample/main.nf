process SEQKIT_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = "fastq"
    if ("${fastx}" ==~ /.+\.(fasta|fa|fas|fna|fsa)(\.gz)?/) {
        extension = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if ("${prefix}.${extension}" == "${fastx}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    seqkit \\
        sample \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fastx} \\
        -o ${prefix}.${extension}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta|.+\.fasta\.gz|.+\.fa|.+\.fa\.gz|.+\.fas|.+\.fas\.gz|.+\.fna|.+\.fna\.gz|.+\.fsa|.+\.fsa\.gz/) {
        extension = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if ("${prefix}.${extension}" == "${fastx}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    def create_cmd = extension.endsWith('.gz') ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${extension}
    """
}
