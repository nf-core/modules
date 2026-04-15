process GENOMAD_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9c/9ce142cdc455bfd9d969463e057da9ee362f7274e6c9fbeb0381c0e3234cae89/data'
        : 'community.wave.seqera.io/library/genomad:1.12.0--17634a7f0b465d30'}"

    output:
    path "genomad_db/" , emit: genomad_db
    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/^.*geNomad, version //; s/ .*//'"), topic: versions, emit: versions_genomad

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    genomad \\
        download-database \\
        ${args} \\
        .
    """

    stub:
    """
    mkdir genomad_db
    touch genomad_db/genomad_db
    touch genomad_db/genomad_db.dbtype
    touch genomad_db/genomad_db.index
    touch genomad_db/genomad_db.lookup
    touch genomad_db/genomad_db.source
    touch genomad_db/genomad_db_h
    touch genomad_db/genomad_db_h.dbtype
    touch genomad_db/genomad_db_h.index
    touch genomad_db/genomad_db_mapping
    touch genomad_db/genomad_db_taxonomy
    touch genomad_db/genomad_integrase_db
    touch genomad_db/genomad_integrase_db.dbtype
    touch genomad_db/genomad_integrase_db.index
    touch genomad_db/genomad_integrase_db.lookup
    touch genomad_db/genomad_integrase_db.source
    touch genomad_db/genomad_integrase_db_h
    touch genomad_db/genomad_integrase_db_h.dbtype
    touch genomad_db/genomad_integrase_db_h.index
    touch genomad_db/genomad_marker_metadata.tsv
    touch genomad_db/genomad_mini_db
    touch genomad_db/genomad_mini_db.dbtype
    touch genomad_db/genomad_mini_db.index
    touch genomad_db/genomad_mini_db.lookup
    touch genomad_db/genomad_mini_db.source
    touch genomad_db/genomad_mini_db_h
    touch genomad_db/genomad_mini_db_h.dbtype
    touch genomad_db/genomad_mini_db_h.index
    touch genomad_db/genomad_mini_db_mapping
    touch genomad_db/genomad_mini_db_taxonomy
    touch genomad_db/mini_set_ids
    touch genomad_db/names.dmp
    touch genomad_db/nodes.dmp
    touch genomad_db/plasmid_hallmark_annotation.txt
    touch genomad_db/version.txt
    touch genomad_db/virus_hallmark_annotation.txt
    """
}
