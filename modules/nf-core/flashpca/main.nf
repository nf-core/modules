process FLASHPCA {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/dbaku42/flashpca:2.0"
    containerOptions "--entrypoint ''"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.pcs.txt")           , emit: pcs,          optional: true
    tuple val(meta), path("*.eigenvectors.txt")  , emit: eigenvectors, optional: true
    tuple val(meta), path("*.eigenvalues.txt")   , emit: eigenvalues,  optional: true
    tuple val(meta), path("*.pve.txt")           , emit: pve,          optional: true
    tuple val(meta), path("*.loadings.txt")      , emit: loadings,     optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L ${bed}  input.bed
    cp -L ${bim}  input.bim
    cp -L ${fam}  input.fam

    flashpca \\
        --bfile input \\
        --numthreads ${task.cpus} \\
        --outpc ${prefix}.pcs.txt \\
        --outvec ${prefix}.eigenvectors.txt \\
        --outval ${prefix}.eigenvalues.txt \\
        --outpve ${prefix}.pve.txt \\
        --outload ${prefix}.loadings.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flashpca: 2.0
    END_VERSIONS
    """
}