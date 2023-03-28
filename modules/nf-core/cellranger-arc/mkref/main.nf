process CELLRANGER_ARC_MKREF {
    tag "$reference_config"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_MKREF module does not support Conda.
        Please use docker or singularity containers."
    }
    container "heylf/cellranger-arc:2.0.2"

    input:
    path fasta
    path gtf
    path motifs
    path reference_config
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellranger-arc \\
        mkref \\
        --config=$reference_config

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
