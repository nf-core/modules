// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ABRICATE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b75992c3e2d82899245f0e5839007cc1dc053cfa72d6018a66729e2ba540af75/data'
        : 'community.wave.seqera.io/library/abricate:1.0.1--cb0de91fb575483f'}"

    input:
    tuple val(meta), path(assembly)
    path databasedir

    output:
    tuple val(meta), path("*.txt"), emit: report
    tuple val("${task.process}"), val('abricate'), eval("abricate --version | sed 's/^.* //' "), emit: versions_abricate, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def datadir = databasedir ? "--datadir ${databasedir}" : ''
    """
    ## Symlink when necessary to rename the file to allow specifying the prefix variable inside report
    ## As the variable is what is used as the sample ID in the report file
    if [[ "${assembly}" != "${prefix}.fasta" ]]; then
        ln -s ${assembly} ${prefix}.fasta
    fi

    abricate \\
        ${prefix}.fasta \\
        ${args} \\
        ${datadir} \\
        --threads ${task.cpus} \\
        > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
