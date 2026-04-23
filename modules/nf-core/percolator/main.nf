process PERCOLATOR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/percolator:3.7.1--h6351f2a_0'
        : 'quay.io/biocontainers/percolator:3.7.1--h6351f2a_0'}"

    input:
    tuple val(meta), path(peptide_identification)

    output:
    tuple val(meta), path("${prefix}.pout.xml"), emit: pout_xml
    tuple val(meta), path("${prefix}.pep.xml"), emit: pout_pepxml
    tuple val(meta), path("${prefix}.features.pin"), emit: features_pin
    tuple val(meta), path("${prefix}.weights.tsv"), emit: weights
    tuple val(meta), path("${prefix}.pep.target.pout"), emit: target_peptides, optional: true
    tuple val(meta), path("${prefix}.pep.decoy.pout"), emit: decoy_peptides, optional: true
    tuple val(meta), path("${prefix}.psm.target.pout"), emit: target_psms
    tuple val(meta), path("${prefix}.psm.decoy.pout"), emit: decoy_psms
    tuple val(meta), path("${prefix}.protein.target.pout"), emit: target_proteins, optional: true
    tuple val(meta), path("${prefix}.protein.decoy.pout"), emit: decoy_proteins, optional: true
    tuple val("${task.process}"), val('percolator'), eval('percolator --help 2>&1 | head -1 | sed "s;Percolator version \\([^,]*\\),.*;\\1;"'), topic: versions, emit: versions_percolator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    percolator \\
        ${args} \\
        --num-threads ${task.cpus} \\
        --xmloutput ${prefix}.pout.xml \\
        --pepxml-output ${prefix}.pep.xml \\
        --tab-out ${prefix}.features.pin \\
        --weights ${prefix}.weights.tsv \\
        --results-peptides ${prefix}.pep.target.pout \\
        --decoy-results-peptides ${prefix}.pep.decoy.pout \\
        --results-psms ${prefix}.psm.target.pout \\
        --decoy-results-psms ${prefix}.psm.decoy.pout \\
        --results-proteins ${prefix}.protein.target.pout \\
        --decoy-results-proteins ${prefix}.protein.decoy.pout \\
        ${peptide_identification}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pout.xml
    touch ${prefix}.pep.xml
    touch ${prefix}.features.pin
    touch ${prefix}.weights.tsv
    touch ${prefix}.pep.target.pout
    touch ${prefix}.pep.decoy.pout
    touch ${prefix}.psm.target.pout
    touch ${prefix}.psm.decoy.pout
    touch ${prefix}.protein.target.pout
    touch ${prefix}.protein.decoy.pout
    """
}
