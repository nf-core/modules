process EXPANSIONHUNTERDENOVO_PROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunterdenovo:0.9.0--hdc99072_3':
        'quay.io/biocontainers/expansionhunterdenovo:0.9.0--hdc99072_3' }"

    input:
    tuple val(meta), path(alignment_file), path(alignment_index)
    tuple val(meta2), path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("*.locus.tsv")        , emit: locus_tsv
    tuple val(meta), path("*.motif.tsv")        , emit: motif_tsv
    tuple val(meta), path("*.str_profile.json") , emit: str_profile
    tuple val("${task.process}"), val('expansionhunterdenovo'), eval("ExpansionHunterDenovo --help |& sed '1!d;s/ExpansionHunter Denovo v//'"), emit: versions_expansionhunterdenovo, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ExpansionHunterDenovo profile \\
        --reads ${alignment_file} \\
        --reference ${fasta} \\
        --output-prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.locus.tsv
    touch ${prefix}.motif.tsv
    touch ${prefix}.str_profile.json
    """
}
