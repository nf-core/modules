process CCNV_REFERENCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ccnv:1.0.0':
        'quay.io/biocontainers/ccnv:1.0.0' }"

    input:
    tuple val(meta), path(cohort_yaml), path(cohort_ccnv, stageAs: 'samples/*'), path(fasta), path(fai)
    tuple val(meta2), path(mappability)
    tuple val(meta3), path(gc_content), path(gc_content_index)

    output:
    tuple val(meta), path("*.ccnv"), emit: ccnv
    tuple val("${task.process}"), val('ccnv'), eval("ccnv --version"), topic: versions, emit: versions_ccnv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gc_cmd = gc_content ? "--gc $gc_content" : ""
    def input_cmd = cohort_yaml && cohort_ccnv ? "--samplesheet $cohort_yaml" : "--fasta ${fasta}"
    def mappability_cmd = mappability ? "--mappability $mappability" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ccnv reference \\
        $args \\
        $gc_cmd \\
        $input_cmd \\
        $mappability_cmd \\
        --prefix $prefix \\
        --threads $task.cpus
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.ccnv
    touch ${prefix}_qc.json
    """
}
