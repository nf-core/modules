process GTDBTK_GTDBTONCBIMAJORITYVOTE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa734cc7e63b0f7d0c04788ec61de5e6a101a07e966d3dde24384d54a9d75e85/data'
        : 'community.wave.seqera.io/library/gtdbtk:2.7.2--64b0fd171db01270'}"

    input:
    tuple val(meta) , path(gtdbtk_outdir), val(gtdbtk_prefix)
    tuple val(meta2), path(ar53_metadata)
    tuple val(meta3), path(bac120_metadata)

    output:
    tuple val(meta), path("*.ncbi.tsv"), emit: tsv
    tuple val("${task.process}"), val('gtdbtk'), eval("gtdbtk --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1") , topic: versions, emit: versions_gtdbtk
    tuple val("${task.process}"), val('gtdb_to_ncbi_majority_vote.py'), eval("gtdb_to_ncbi_majority_vote.py -v 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -n 1"), topic: versions, emit: versions_gtdbtoncbimajorityvote

    when:
    task.ext.when == null || task.ext.when

    script:
    if(!bac120_metadata && !ar53_metadata) {
        log.error("ERROR: Neither the ar53 or bac120 metadata files were provided to GTDBTK_GTDBTONCBIMAJORITYVOTE!")
    }
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def args          = task.ext.args   ?: ""
    def prefix_arg    = gtdbtk_prefix   ? "--gtdbtk_prefix ${gtdbtk_prefix}"          : ""
    bac120            = bac120_metadata ? "--bac120_metadata_file ${bac120_metadata}" : ""
    ar53              = ar53_metadata   ? "--ar53_metadata_file ${ar53_metadata}"     : ""
    """
    gtdb_to_ncbi_majority_vote.py \\
        --gtdbtk_output_dir ${gtdbtk_outdir} \\
        ${prefix_arg} \\
        ${bac120} \\
        ${ar53} \\
        --output_file ${prefix}.ncbi.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gtdbtk.${prefix}.ncbi.tsv
    """
}
