process LEARNMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/learnmsa:2.0.16--pyhdfd78af_0'
        : 'quay.io/biocontainers/learnmsa:2.0.16--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    tuple val("${task.process}"), val('learnmsa'), eval("learnMSA -h | sed -nE 's/.*version ([0-9.]+).*/\\1/p'"), emit: versions_learnmsa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("LearnMSA align module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    """
    learnMSA \\
        -i ${fasta} \\
        -o "${prefix}.aln" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln
    """
}
