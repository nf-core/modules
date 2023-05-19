
process PLINK_INDEP {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed),  path(bim), path(fam)
    val(window_size)
    val(variant_count)
    val(variance_inflation_factor)

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
        --threads $task.cpus \\
        --indep ${window_size} ${variant_count} ${variance_inflation_factor} \\
        $args \\
        --out $prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
