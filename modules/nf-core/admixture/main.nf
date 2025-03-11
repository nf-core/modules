process ADMIXTURE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/admixture:1.3.0--0':
        'biocontainers/admixture:1.3.0--0' }"

    input:
    tuple val(meta), path (bed_ped_geno), path(bim_map), path(fam)
    val K


    output:
    tuple val(meta), path("*.Q")    , emit: ancestry_fractions
    tuple val(meta), path("*.P")    , emit: allele_frequencies
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    admixture \\
        $bed_ped_geno \\
        $K \\
        -j$task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(echo \$(admixture 2>&1) | head -n 1 | grep -o "ADMIXTURE Version [0-9.]*" | sed 's/ADMIXTURE Version //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.Q"
    touch "${prefix}.P"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(echo \$(admixture 2>&1) | head -n 1 | grep -o "ADMIXTURE Version [0-9.]*" | sed 's/ADMIXTURE Version //' )
    END_VERSIONS
    """
}
