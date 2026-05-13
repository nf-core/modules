process DEEPMASED_FEATURES {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deepmased:0.3.1--pyh5ca1d4c_0':
        'quay.io/biocontainers/deepmased:0.3.1--pyh5ca1d4c_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("${prefix}_feature_file_paths.tsv"), path("*_feats.tsv"), emit: features
    tuple val("${task.process}"), val('deepmased'), val('0.3.1'), topic: versions, emit: versions_deepmased_features
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "bam\\tfasta" > ${prefix}_file_paths.tsv
    echo -e "${bam}\\t${fasta}" >> ${prefix}_file_paths.tsv

    DeepMAsED features \\
        ${prefix}_file_paths.tsv \\
        -p ${task.cpus} \\
        -o . \\
        -n ${prefix}_feature_file_paths.tsv \\
        ${args}
    """

    stub:
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_feature_file_paths.tsv
    touch ${prefix}_feats.tsv
    """
}
