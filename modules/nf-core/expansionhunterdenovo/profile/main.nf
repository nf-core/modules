process EXPANSIONHUNTERDENOVO_PROFILE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::expansionhunterdenovo=0.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunterdenovo:0.9.0--hdc99072_3':
        'biocontainers/expansionhunterdenovo:0.9.0--hdc99072_3' }"

    input:
    tuple val(meta), path(alignment_file), path(alignment_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.locus.tsv")        , emit: locus_tsv
    tuple val(meta), path("*.motif.tsv")        , emit: motif_tsv
    tuple val(meta), path("*.str_profile.json") , emit: str_profile
    path "versions.yml"                         , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunterdenovo: \$(echo \$(ExpansionHunterDenovo --help 2>&1) | sed -e "s/ExpansionHunter Denovo v//;s/ Usage.*//")
    END_VERSIONS
    """
}
