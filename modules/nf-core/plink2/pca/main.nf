process PLINK2_PCA {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a2.3--h712d239_1' :
        'biocontainers/plink2:2.00a2.3--h712d239_1' }"

    input:
    tuple val(meta), val(npcs), val(use_approx), path(pgen), path(psam), path(pvar)

    output:
    tuple val(meta), path("*.eigenvec")    , emit: eigenvec
    tuple val(meta), path("*.eigenval")    , emit: eigenval
    tuple val(meta), path("*.log")         , emit: logfile
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def approx_option = use_approx ? "approx" : ""
    def n_pcs = npcs ? npcs : 10
    """
    plink2 \\
        --pca ${n_pcs} ${approx_option} \\
        --memory ${task.memory.toMega()} \\
        $args \\
        --threads $task.cpus \\
        --pfile ${pgen.baseName} \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plinkpca: \$(plink2 --version |& sed '1!d ; s/plink2 //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eigenvec ${prefix}.eigenval

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
