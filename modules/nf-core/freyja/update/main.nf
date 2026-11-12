// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FREYJA_UPDATE {
    tag "$db_name"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f0b023f9c9e25f338d72ed6710e4d2a54a4c3a9a3ac34f17e08683894dc40830/data'
:         'community.wave.seqera.io/library/freyja:2.0.6--be3cbd6f14f9d542' }"

    input:
    val db_name

    output:
    path "${db_name}/*barcodes.*"            , emit: barcodes
    path "${db_name}/*lineages.yml"          , emit: lineages_topology
    path "${db_name}/*pathogen_config.yml"   , emit: config
    path "${db_name}/*curated_lineages.json" , emit: lineages_meta, optional: true
    tuple val("${task.process}"), val('freyja'), eval("freyja --version | sed 's/.* //'"), topic: versions, emit: versions_freyja

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p $db_name
    freyja \\
        update \\
        $args \\
        --outdir $db_name
    """

    stub:
    """
    mkdir $db_name

    touch "${db_name}/usher_barcodes.csv"
    touch "${db_name}/lineages.yml"
    touch "${db_name}/curated_lineages.json"
    touch "${db_name}/pathogen_config.yml"

    """
}
