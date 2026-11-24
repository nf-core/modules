process CCNV_REFERENCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cd/cd552e195217d37801cf86348dbe875fa66ae2014147e481fd6be6e80dbbcd90/data'
        : 'community.wave.seqera.io/library/ccnv:1.1.0--726d8b142feb0e3d'}"

    input:
    tuple val(meta), path(cohort_yaml), path(cohort_ccnv), path(fasta), path(fai)
    tuple val(meta2), path(mappability)
    tuple val(meta3), path(gc_content), path(gc_content_index)

    output:
    tuple val(meta), path("*.ccnv"), emit: ccnv
    tuple val(meta), path("*_qc.json"), emit: qc_json, optional: true
    tuple val("${task.process}"), val('ccnv'), eval("ccnv --version 2>&1 | sed 's/ccnv //'"), topic: versions, emit: versions_ccnv

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
