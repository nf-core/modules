process PURGEDUPS_GETSEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--py39h7132678_1':
        'biocontainers/purge_dups:1.2.6--py39h7132678_1' }"

    input:
    tuple val(meta), path(assembly), path(bed)

    output:
    tuple val(meta), path("*.hap.fa")   , emit: haplotigs
    tuple val(meta), path("*.purged.fa"), emit: purged
    // WARN: Incorrect version printed inside the container, please check this if bumping version ( \$( purge_dups -h |& sed '3!d; s/.*: //' ))
    tuple val("${task.process}"), val('purge_dups'), val('1.2.6'), emit: versions_purgedups, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // By default, we include the '-e' option in task.ext.args as this is considered the 'safe default'
    // Check if we have a CharSequence (both String and GString inherit) and thus we can apply an empty string, which is otherwise 'falsy'
    // See https://github.com/dfguan/purge_dups?tab=readme-ov-file#step-3-get-purged-primary-and-haplotig-sequences-from-draft-assembly
    def args = (task.ext.args instanceof CharSequence) ? task.ext.args : '-e'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_seqs \\
        ${args} \\
        ${bed} \\
        -p ${prefix} \\
        ${assembly}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.hap.fa"
    touch "${prefix}.purged.fa"
    """
}
