process EXPANSIONHUNTERDENOVO_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunterdenovo:0.9.0--hdc99072_3':
        'quay.io/biocontainers/expansionhunterdenovo:0.9.0--hdc99072_3' }"

    input:
    tuple val(meta), path(manifest)
    tuple val(meta2), path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("*.multisample_profile.json"), emit: merged_profiles
    tuple val("${task.process}"), val('expansionhunterdenovo'), eval("ExpansionHunterDenovo --help |& sed '1!d;s/ExpansionHunter Denovo v//'"), emit: versions_expansionhunterdenovo, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.multisample_profile.json
    """
}
