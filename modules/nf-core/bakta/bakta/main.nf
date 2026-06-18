process BAKTA_BAKTA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50b75335f6394ae83fd05f364db27ee2eb75f4170e3525bb2aea47ad717a9e64/data'
        : 'community.wave.seqera.io/library/bakta_diamond:7830b94718da4f96'}"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path regions
    path hmms

    output:
    tuple val(meta), path("${prefix}.embl"), emit: embl
    tuple val(meta), path("${prefix}.faa"), emit: faa
    tuple val(meta), path("${prefix}.ffn"), emit: ffn
    tuple val(meta), path("${prefix}.fna"), emit: fna
    tuple val(meta), path("${prefix}.gbff"), emit: gbff
    tuple val(meta), path("${prefix}.gff3"), emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('bakta'), eval("bakta --version 2>&1 | sed 's/.*bakta //'"), emit: versions_bakta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    def regions_opt = regions ? "--regions ${regions}" : ""
    def hmms_opt = hmms ? "--hmms ${hmms}" : ""

    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    bakta \\
        ${fasta} \\
        ${args} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${proteins_opt} \\
        ${prodigal_tf_opt} \\
        ${regions_opt} \\
        ${hmms_opt} \\
        --db ${db}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    touch ${prefix}.embl
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.fna
    touch ${prefix}.gbff
    touch ${prefix}.gff3
    touch ${prefix}.hypotheticals.tsv
    touch ${prefix}.hypotheticals.faa
    touch ${prefix}.tsv
    touch ${prefix}.txt
    """
}
