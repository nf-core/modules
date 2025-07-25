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
    tuple val(meta), path("${prefix}.vamb/bins/*.fna.gz")             , emit: bins             , optional: true
    tuple val(meta), path("${prefix}.vamb/vae*_clusters_metadata.tsv"), emit: clusters_metadata
    tuple val(meta), path("${prefix}.vamb/vae*_clusters_split.tsv")   , emit: clusters_split   , optional: true
    tuple val(meta), path("${prefix}.vamb/vae*_clusters_unsplit.tsv") , emit: clusters_unsplit
    tuple val(meta), path("${prefix}.vamb/results_taxometer.tsv")     , emit: taxometer_results, optional: true
    tuple val(meta), path("${prefix}.vamb/latent.npz")                , emit: latent_encoding  , optional: true
    tuple val(meta), path("${prefix}.vamb/abundance.npz")             , emit: abundance
    tuple val(meta), path("${prefix}.vamb/composition.npz")           , emit: composition
    tuple val(meta), path("${prefix}.vamb/log.txt")                   , emit: log
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if(bams && abundance_tsv) {
        error("ERROR: Both bams and abundance TSV supplied to Vamb!")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode = taxonomy ? "taxvamb" : "default"
    def depth_input = abundance_tsv ? "--abundance_tsv ${abundance_tsv}" : "--bamdir bams/"
    def tax_input = taxonomy ? "--taxonomy ${taxonomy}" : ""
    """
    vamb bin \\
        ${mode} \\
        -p ${task.cpus} \\
        --outdir ${prefix}.vamb/ \\
        --fasta ${assembly} \\
        ${depth_input} \\
        ${tax_input} \\
        ${args}

    find ${prefix}.vamb/bins -name "*.fna" -exec gzip {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """

    stub:
    if(bams && abundance_tsv) {
        error("ERROR: Both bams and abundance TSV supplied to Vamb!")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.vamb/bins

    echo "" | gzip > ${prefix}.vamb/bins/1.fna.gz
    echo "" | gzip > ${prefix}.vamb/bins/2.fna.gz

    touch ${prefix}.vamb/results_taxometer.tsv
    touch ${prefix}.vamb/predictor_model.pt
    touch ${prefix}.vamb/vae_clusters_metadata.tsv
    touch ${prefix}.vamb/vae_clusters_split.tsv
    touch ${prefix}.vamb/vae_clusters_unsplit.tsv
    touch ${prefix}.vamb/latent.npz
    touch ${prefix}.vamb/model.pt
    touch ${prefix}.vamb/abundance.npz
    touch ${prefix}.vamb/composition.npz
    touch ${prefix}.vamb/log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """
}
