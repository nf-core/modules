process RIPGREP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/deecf2e0aab3f4ba7c04eb24724f5c80550e39bec80760b18c7b52f8efa8e6b3/data':
        'community.wave.seqera.io/library/pigz_ripgrep:94e5407412b666ab' }"

    input:
    tuple val(meta), path(files, arity: '1..*')
    val pattern
    path pattern_file
    val compress

    output:
    tuple val(meta), path("*.txt{.gz,}"), emit: txt
    tuple val("${task.process}"), val('ripgrep'), eval("rg --version |& sed '1!d ; s/ripgrep // ; s/ .*//'"),  emit: versions_ripgrep, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.txt.gz" : "> ${prefix}.txt"
    if (pattern && pattern_file)    {
        error("RIPGREP provided both a pattern and a pattern file!")
    }
    if (!compress && files.contains("${prefix}.txt") || compress && files.contains("${prefix}.txt.gz")) {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    def pattern_arg = pattern_file ? "-f ${pattern_file}" : "${pattern}"
    """
    rg \\
        $args \\
        --threads $task.cpus \\
        ${pattern_arg} \\
        ${files} \\
        ${write_output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_output = compress ? "echo | gzip > ${prefix}.txt.gz" : "touch ${prefix}.txt"
    """
    $write_output
    """
}
