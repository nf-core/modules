process GTDBTK_GTDBTONCBIMAJORITYVOTE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.4.1--pyhdfd78af_1':
        'biocontainers/gtdbtk:2.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta) , path(gtdbtk_output)
    tuple val(meta2), path(ar53_metadata)
    tuple val(meta3), path(bac120_metadata)

    output:
    tuple val(meta), path("*.ncbi.tsv"), emit: tsv
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gtdbtk_prefix = gtdbtk_output.find { it.name =~ /summary\.tsv/}.getName() - ~/\.(bac120|ar53)\.summary\.tsv/
    """
    gtdb_to_ncbi_majority_vote.py \\
        --gtdbtk_output_dir . \\
        ${bac120_metadata} \\
        ${ar53_metadata} \\
        --gtdbtk_prefix ${gtdbtk_prefix} \\
        --output_file gtdbtk.${prefix}.ncbi.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(gtdb_to_ncbi_majority_vote.py -v 2>/dev/null | head -n 1 | grep -Po "(\\d+\\.)+\\d+")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gtdbtk.${prefix}.ncbi.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(gtdb_to_ncbi_majority_vote.py -v 2>/dev/null | head -n 1 | grep -Po "(\\d+\\.)+\\d+")
    END_VERSIONS
    """
}
