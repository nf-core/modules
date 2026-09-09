process PCAONE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f908496d8e5dc44b1f5caa8d330e714250aacc9f65c1577791436914897af56/data':
        'community.wave.seqera.io/library/pcaone:0.7.2--22489ef995a570f0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.eigvecs")   , emit: eigvecs
    tuple val(meta), path("*.eigvecs2")  , emit: eigvecs2
    tuple val(meta), path("*.eigvals")   , emit: eigvals
    tuple val(meta), path("*.loadings")  , emit: loadings  , optional: true
    tuple val(meta), path("*.mbim")      , emit: mbim      , optional: true
    tuple val(meta), path("*.log")       , emit: log
    tuple val("${task.process}"), val('pcaone'), eval('PCAone 2>&1 | grep -m1 -oE "v[0-9]+\\.[0-9]+\\.[0-9]+" | tr -d v'), topic: versions, emit: versions_pcaone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // bed/bim/fam must share the same file stem for PCAone's --bfile option
    def bfile  = bed.baseName
    """
    PCAone \\
        --bfile ${bfile} \\
        --threads ${task.cpus} \\
        --out ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eigvecs
    touch ${prefix}.eigvecs2
    touch ${prefix}.eigvals
    touch ${prefix}.log
    """
}
