process FASTAGUARD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastaguard:0.6.0--hfa8f182_0':
        'quay.io/biocontainers/fastaguard:0.6.0--hfa8f182_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fastaguard.html"), emit: html
    tuple val(meta), path("*.fastaguard.json"), emit: json
    tuple val(meta), path("*.fastaguard.tsv"), emit: tsv
    tuple val(meta), path("*.fastaguard_mqc.json"), emit: mqc
    tuple val("${task.process}"), val('fastaguard'), eval('fastaguard --version | cut -d " " -f 2'), emit: versions_fastaguard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    fastaguard ${fasta} \
      ${args} \
      --out ${prefix}.fastaguard.html \
      --json ${prefix}.fastaguard.json \
      --tsv ${prefix}.fastaguard.tsv \
      --multiqc ${prefix}.fastaguard_mqc.json
    """
}
