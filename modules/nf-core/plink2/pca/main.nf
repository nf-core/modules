process PLINK2_PCA {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a2.3--h712d239_1' :
        'biocontainers/plink2:2.00a2.3--h712d239_1' }"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)
    val(bed)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple path("*.eigenvec"), path("*.eigenval"), path("*.log"), emit: eigen
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink2 \\
        --pca \\
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
    touch ${prefix}.eigenvec ${prefix}.eigenval ${prefix}.log 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plinkpca: \$(plink2 --version |& sed '1!d ; s/plink2 //')
    END_VERSIONS
    """
}
