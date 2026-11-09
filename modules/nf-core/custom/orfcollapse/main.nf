// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_ORFCOLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/7259cf5c00cd90b521cdda013955b69043542ad2b93f786c6ed53ecfc04a8422/data'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml:26e261815660ec9d'}"

    input:
    tuple val(meta), path(bed12, stageAs: 'input/*'), path(catalogue_tsv, stageAs: 'input/*'), path(orf_to_gene_tsv, stageAs: 'input/*'), path(aa_fasta, stageAs: 'input/*'), path(cluster_tsv, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.bed12")                     , emit: bed12
    tuple val(meta), path("${prefix}.tsv")                       , emit: catalogue_tsv
    tuple val(meta), path("${prefix}.orf_to_gene.tsv")           , emit: orf_to_gene_tsv
    tuple val(meta), path("${prefix}.mqc.tsv")                   , emit: multiqc
    tuple val(meta), path("${prefix}.fasta")                     , emit: aa_fasta
    tuple val(meta), path("${prefix}.consensus.bed12")           , emit: consensus_bed12
    tuple val(meta), path("${prefix}.consensus.tsv")             , emit: consensus_tsv
    tuple val(meta), path("${prefix}.consensus.orf_to_gene.tsv") , emit: consensus_orf_to_gene_tsv
    path "versions.yml"                                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.catalogue"
    args   = task.ext.args ?: ''

    """
    echo ${args}
    """

    template 'orfcollapse.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.catalogue"
    """
    touch ${prefix}.bed12 ${prefix}.tsv ${prefix}.orf_to_gene.tsv ${prefix}.mqc.tsv ${prefix}.fasta
    touch ${prefix}.consensus.bed12 ${prefix}.consensus.tsv ${prefix}.consensus.orf_to_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """
}
