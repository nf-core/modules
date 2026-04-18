process VAMB_BIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vamb:5.0.4--pyhdfd78af_0':
        'biocontainers/vamb:5.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(abundance_tsv), path(bams, stageAs: "bams/*"), path(taxonomy)

    output:
    tuple val(meta), path("${prefix}/bins/*.fna.gz")             , emit: bins             , optional: true
    tuple val(meta), path("${prefix}/vae*_clusters_metadata.tsv"), emit: clusters_metadata
    tuple val(meta), path("${prefix}/vae*_clusters_split.tsv")   , emit: clusters_split   , optional: true
    tuple val(meta), path("${prefix}/vae*_clusters_unsplit.tsv") , emit: clusters_unsplit
    tuple val(meta), path("${prefix}/results_taxometer.tsv")     , emit: taxometer_results, optional: true
    tuple val(meta), path("${prefix}/latent.npz")                , emit: latent_encoding  , optional: true
    tuple val(meta), path("${prefix}/abundance.npz")             , emit: abundance
    tuple val(meta), path("${prefix}/composition.npz")           , emit: composition
    tuple val(meta), path("${prefix}/log.txt")                   , emit: log
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if(bams && abundance_tsv) {
        error("ERROR: Both bams and abundance TSV supplied to Vamb! Please only supply one.")
    }
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def mode    = taxonomy ? "taxvamb" : "default"
    depth_input = abundance_tsv ? "--abundance_tsv ${abundance_tsv}" : "--bamdir bams/"
    tax_input   = taxonomy ? "--taxonomy ${taxonomy}" : ""
    """
    vamb bin \\
        ${mode} \\
        -p ${task.cpus} \\
        --outdir ${prefix}/ \\
        --fasta ${assembly} \\
        ${depth_input} \\
        ${tax_input} \\
        ${args}

    find ${prefix}/bins -name "*.fna" -exec gzip {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """

    stub:
    if(bams && abundance_tsv) {
        error("ERROR: Both bams and abundance TSV supplied to Vamb! Please only supply one.")
    }
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/bins

    echo "" | gzip > ${prefix}/bins/1.fna.gz
    echo "" | gzip > ${prefix}/bins/2.fna.gz

    touch ${prefix}/results_taxometer.tsv
    touch ${prefix}/predictor_model.pt
    touch ${prefix}/vae_clusters_metadata.tsv
    touch ${prefix}/vae_clusters_split.tsv
    touch ${prefix}/vae_clusters_unsplit.tsv
    touch ${prefix}/latent.npz
    touch ${prefix}/model.pt
    touch ${prefix}/abundance.npz
    touch ${prefix}/composition.npz
    touch ${prefix}/log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """
}
