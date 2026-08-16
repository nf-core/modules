process FOLDSEEK_EASYSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
        'quay.io/biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"

    input:
    tuple val(meta) , path(pdb)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("${prefix}.m8"), emit: aln
    tuple val("${task.process}"), val('foldseek'), eval("foldseek --help |& sed -n 's/.*Version: //p'"), emit: versions_foldseek, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    """
    foldseek \\
        easy-search \\
        ${pdb} \\
        ${db}/${prefix2} \\
        ${prefix}.m8 \\
        tmpFolder \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.m8
    """
}
