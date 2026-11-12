process BULLSEYE_FINDEDITSITES {
    tag "$meta.id vs $control_meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d1b3ef00779596d2255fc5dffeae072301ccc37d6cdb941c2aa8d67534ced41/data'
        : 'community.wave.seqera.io/library/bullseye:1.1.0--6a0a52edc73ac840'}"

    input:
    tuple val(meta), path(edited_matrix), path(edited_tbi), val(control_meta), path(control_matrix), path(control_tbi)
    path annotation
    val edit_type
    val min_edit
    val max_edit
    val edit_fold_threshold
    val min_edit_sites
    val edited_min_coverage
    val control_min_coverage

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('bullseye'), eval("grep -m1 -i 'version' \$(command -v Find_edit_sites.pl) | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?'"), topic: versions, emit: versions_bullseye

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.vs.${control_meta.id}"
    """
    perl -S Find_edit_sites.pl \\
        --annotationFile ${annotation} \\
        --EditedMatrix ${edited_matrix} \\
        --controlMatrix ${control_matrix} \\
        --editType ${edit_type} \\
        --minEdit ${min_edit} \\
        --maxEdit ${max_edit} \\
        --editFoldThreshold ${edit_fold_threshold} \\
        --MinEditSites ${min_edit_sites} \\
        --EditedMinCoverage ${edited_min_coverage} \\
        --ControlMinCoverage ${control_min_coverage} \\
        --outfile ${prefix}.bed \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.vs.${control_meta.id}"
    """
    touch ${prefix}.bed
    """
}
