process DOTMATCH_CRISPR_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dotmatch:0.2.1--py314h118bc1c_0' :
        'quay.io/biocontainers/dotmatch:0.2.1--py314h118bc1c_0' }"

    input:
    tuple val(meta), path(reads)
    path library
    val guide_start
    val guide_length
    val k
    val metric

    output:
    tuple val(meta), path("*.counts.mageck.tsv"), emit: counts
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.sample_qc.tsv"), emit: sample_qc
    tuple val("${task.process}"), val('dotmatch'), eval('dotmatch --version | sed "s/^dotmatch //"'), emit: versions_dotmatch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample_id\\tfastq\\n' > samples.tsv
    printf '%s\\t%s\\n' '${meta.id}' '${reads}' >> samples.tsv

    dotmatch crispr-count \\
        --library ${library} \\
        --samples samples.tsv \\
        --guide-start ${guide_start} \\
        --guide-length ${guide_length} \\
        --k ${k} \\
        --metric ${metric} \\
        --ambiguity-policy radius \\
        --threads ${task.cpus} \\
        --out ${prefix}.counts.mageck.tsv \\
        --summary ${prefix}.summary.json \\
        --sample-qc ${prefix}.sample_qc.tsv \\
        --ambiguous discard \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.mageck.tsv
    touch ${prefix}.summary.json
    touch ${prefix}.sample_qc.tsv
    """
}
