process LTRRETRIEVER_LAI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_retriever:3.0.5--hdfd78af_0':
        'quay.io/biocontainers/ltr_retriever:3.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path pass_list
    path annotation_out
    path monoploid_seqs

    output:
    tuple val(meta), path("*.LAI.log")  , emit: log
    tuple val(meta), path("*.LAI.out")  , emit: lai_out     , optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val("lai"), val("beta3.2"), emit: versions_lai, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args     ?: ''
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def monoploid_param = monoploid_seqs    ? "-mono ${monoploid_seqs}"                       : ''
    def lai_output_name = monoploid_seqs    ? "${annotation_out}.${monoploid_seqs}.out.LAI" : "${annotation_out}.LAI"
    """
    LAI \\
        -genome ${fasta} \\
        -intact ${pass_list} \\
        -all ${annotation_out} \\
        -t ${task.cpus} \\
        ${monoploid_param} \\
        ${args} \\
        >| >(tee "${prefix}.LAI.log") \\
        || echo "LAI failed! See ${prefix}.LAI.log"

    mv \\
        ${lai_output_name} \\
        "${prefix}.LAI.out" \\
        || echo "LAI failed to estimate assembly index. See ${prefix}.LAI.log"
    """

    stub:
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    """
    touch "${prefix}.LAI.log"
    touch "${prefix}.LAI.out"
    """
}
