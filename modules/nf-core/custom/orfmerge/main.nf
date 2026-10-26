process CUSTOM_ORFMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/caf07471ede03189e710de9163ae02e17bf33643bdcd6f1084ef57484174bf7c/data'
        : 'community.wave.seqera.io/library/python_pyyaml:14e95bf7fff3fa6c'}"

    input:
    tuple val(meta), path(bed12s, arity: '1..*', stageAs: 'beds/*'), path(tsvs, arity: '1..*', stageAs: 'tsvs/*')

    output:
    tuple val(meta), path("${prefix}.bed12")                     , emit: bed12
    tuple val(meta), path("${prefix}.tsv")                       , emit: catalogue_tsv
    tuple val(meta), path("${prefix}.orf_to_gene.tsv")           , emit: orf_to_gene_tsv
    tuple val(meta), path("${prefix}.mqc.tsv")                   , emit: multiqc
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

    template 'orfmerge.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.catalogue"
    """
    touch ${prefix}.bed12
    touch ${prefix}.tsv
    touch ${prefix}.orf_to_gene.tsv
    touch ${prefix}.mqc.tsv
    touch ${prefix}.consensus.bed12
    touch ${prefix}.consensus.tsv
    touch ${prefix}.consensus.orf_to_gene.tsv

    python -c "import platform, yaml; yaml.safe_dump({'${task.process}': {'python': platform.python_version()}}, open('versions.yml', 'w'), default_flow_style=False, sort_keys=False)"
    """
}
