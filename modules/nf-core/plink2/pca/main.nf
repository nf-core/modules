process PLINK2_PCA {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dfbbefb7c4dccd8c211a1c8bf94988dcd6a58a147e5c22e2553e80af907569/data'
:         'community.wave.seqera.io/library/plink2:2.00a5.10--223badfc86e2b6e6' }"

    input:
    tuple val(meta), val(npcs), val(use_approx), path(pgen), path(psam), path(pvar)

    output:
    tuple val(meta), path("*.eigenvec"), emit: evecfile
    tuple val(meta), path("*.eigenval"), emit: evfile
    tuple val(meta), path("*.log"), emit: logfile
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def approx_option = use_approx ? "approx" : ""
    def n_pcs = npcs ?: 10
    """
    plink2 \\
        --pca ${n_pcs} ${approx_option} \\
        --memory ${task.memory.toMega()} \\
        $args \\
        --threads $task.cpus \\
        --pfile ${pgen.baseName} \\
        --out ${prefix}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eigenvec
    touch ${prefix}.eigenval
    touch ${prefix}.log
    """
}
