process EXPANSIONHUNTERDENOVO_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::expansionhunterdenovo=0.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunterdenovo:0.9.0--hdc99072_3':
        'quay.io/biocontainers/expansionhunterdenovo:0.9.0--hdc99072_3' }"

    input:
    tuple val(meta), path(manifest)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.multisample_profile.json"), emit: merged_profiles
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ExpansionHunterDenovo merge \\
        --manifest ${manifest} \\
        --reference ${fasta} \\
        --output-prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunterdenovo: \$(echo \$(ExpansionHunterDenovo --help 2>&1) | sed -e "s/ExpansionHunter Denovo v//;s/ Usage.*//")
    END_VERSIONS
    """
}
