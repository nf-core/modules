process GEDI_INDEXGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/bae29fa913dea79a3dcdbfbf544f0391f82bbfdbf3e6430f71db45ba21d6cf79/data' :
        'community.wave.seqera.io/library/gedi_indexgenome:cfca16738f306c86' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path("price_index")                                                  , emit: index
    tuple val("${task.process}"), val('gedi'), eval("gedi -e Version 2>&1 | sed -n 's/.*Gedi version \\([^ ]*\\).*/\\1/p' | head -n 1"), topic: versions, emit: versions_gedi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = meta.id ?: 'reference'
    """
    mkdir -p price_index
    gedi -e IndexGenome \\
        -s ${fasta} \\
        -a ${gtf} \\
        -n ${name} \\
        -f price_index \\
        -o price_index/${name}.oml \\
        -nomapping \\
        -p \\
        ${args}
    """

    stub:
    def name = meta.id ?: 'reference'
    """
    mkdir -p price_index
    touch price_index/${name}.oml
    """
}
