process PLINK2_PCAWTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0'
        : 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'}"

    input:
    tuple val(meta), val(npcs), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.eigenvec.allele"), emit: allele_weights
    tuple val(meta), path("*.eigenvec"), emit: eigenvec
    tuple val(meta), path("*.eigenval"), emit: eigenval
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_pcs = npcs ?: 10
    """
    plink2 \\
        --pgen ${pgen} \\
        --pvar ${pvar} \\
        --psam ${psam} \\
        --pca ${n_pcs} allele-wts \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eigenvec.allele
    touch ${prefix}.eigenvec
    touch ${prefix}.eigenval
    touch ${prefix}.log
    """
}
