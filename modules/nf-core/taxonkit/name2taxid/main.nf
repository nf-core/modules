// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAXONKIT_NAME2TAXID {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91013897d309411f36f0bf958988adf50d8abeeecd910b7debd446bfc1e3f2d8/data'
        : 'community.wave.seqera.io/library/taxonkit:0.15.1--96f7341bd5c149a7'}"

    input:
    tuple val(meta), val(name), path(names_txt)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('taxonkit'), eval("taxonkit version | sed 's/.* v//'"), emit: versions_taxonkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert (!name && names_txt) || (name && !names_txt)
    """
    taxonkit \\
        name2taxid \\
        ${args} \\
        --data-dir ${taxdb} \\
        --threads ${task.cpus} \\
        --out-file ${prefix}.tsv \\
        ${name? "<<< '${name}'": names_txt}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
