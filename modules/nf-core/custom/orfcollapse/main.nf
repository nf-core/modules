process CUSTOM_ORFCOLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1019bd22c111267bcb670fdb128829776f0ca6adfa7b0e2d126f91577d08e3/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:75514f9f977be607' }"

    input:
    tuple val(meta), path(bed12, stageAs: 'input/*'), path(catalogue_tsv, stageAs: 'input/*'), path(orf_to_gene_tsv, stageAs: 'input/*'), path(aa_fasta, stageAs: 'input/*'), path(cluster_tsv, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.bed12")           , emit: bed12
    tuple val(meta), path("${prefix}.tsv")             , emit: catalogue_tsv
    tuple val(meta), path("${prefix}.orf_to_gene.tsv") , emit: orf_to_gene_tsv
    tuple val(meta), path("${prefix}.mqc.tsv")         , emit: multiqc
    tuple val(meta), path("${prefix}.fasta")           , emit: aa_fasta
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.catalogue"
    args   = task.ext.args ?: ''
    template 'orfcollapse.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.catalogue"
    """
    touch ${prefix}.bed12 ${prefix}.tsv ${prefix}.orf_to_gene.tsv ${prefix}.mqc.tsv ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """
}
