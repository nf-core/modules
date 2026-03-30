process PERCOLATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/percolator:3.7.1--h6351f2a_0'
        : 'biocontainers/percolator:3.7.1--h6351f2a_0'}"

    input:
    tuple val(meta), path(peptide_identification)

    output:
    tuple val(meta), path("*.target.pout"), emit: target_pout
    tuple val(meta), path("*.decoy.pout"), emit: decoy_pout

    tuple val("${task.process}"), val('percolator'), eval('percolator --help 2>&1 | head -1 | sed "s;Percolator version \\([^,]*\\),.*;\\1;"'), topic: versions, emit: versions_percolator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    percolator \\
        ${args} \\
        --num-threads ${task.cpus} \\
        --only-psms \\
        --post-processing-tdc \\
        --search-input concatenated \\
        --results-psms ${prefix}.psm.target.pout \\
        --decoy-results-psms ${prefix}.psm.decoy.pout \\
        ${peptide_identification}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.psm.target.pout
    touch ${prefix}.psm.decoy.pout
    """
}
