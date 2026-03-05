process GTDBTK_GTDBTONCBIMAJORITYVOTE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.5.2--pyh1f0d9b5_0':
        'biocontainers/gtdbtk:2.5.2--pyh1f0d9b5_0' }"

    input:
    tuple val(meta) , path(gtdbtk_outdir), val(gtdbtk_prefix)
    tuple val(meta2), path(ar53_metadata)
    tuple val(meta3), path(bac120_metadata)

    output:
    tuple val(meta), path("*.ncbi.tsv"), emit: tsv
    path "versions.yml"                , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py -v 2>/dev/null) | grep -o -E "[0-9]+(\\.[0-9]+)+" | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gtdbtk.${prefix}.ncbi.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py -v 2>/dev/null) | grep -o -E "[0-9]+(\\.[0-9]+)+" | head -n 1)
    END_VERSIONS
    """
}
