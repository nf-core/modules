process PLINK_INDEPPAIRWISE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed),  path(bim), path(fam)
    val(window_size)
    val(variant_count)
    val(r2_threshold)

    output:
    tuple val(meta), path("*.prune.in")                    , emit: prunein
    tuple val(meta), path("*.prune.out")    , optional:true, emit: pruneout
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink \\
        --bed ${bed}  \\
        --bim ${bim}  \\
        --fam ${fam}  \\
        --threads ${task.cpus} \\
        --indep-pairwise ${window_size} ${variant_count} ${r2_threshold} \\
        ${args} \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${meta.id}.prune.in

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
