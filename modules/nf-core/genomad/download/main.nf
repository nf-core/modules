process GENOMAD_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83e31b082e82e714b01040adf90f865fdd229b6d9d926525b821d207271a3922/data'
        : 'community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5'}"

    output:
    path "genomad_db/", emit: genomad_db
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
