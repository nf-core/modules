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
    tuple val(meta), path("${prefix}.pout.xml"), emit: pout_xml, optional: true
    tuple val(meta), path("${prefix}.pep.xml"), emit: pout_pepxml, optional: true
    tuple val(meta), path("${prefix}.features.pin"), emit: features_pin, optional: true
    tuple val(meta), path("${prefix}.weights.tsv"), emit: weights, optional: true
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
    // enhance the helper-flags
    def write_xml = args.contains("--write_xml") ? "--xmloutput ${prefix}.pout.xml" : ""
    def write_pepxml = args.contains("--write_pepxml") ? "--pepxml-output ${prefix}.pep.xml" : ""
    def write_tabout = args.contains("--write_tabout") ? "--tab-out ${prefix}.features.pin" : ""
    def write_weights = args.contains("--write_weights") ? "--weights ${prefix}.weights.tsv" : ""
    def write_peptide_results = args.contains("--write_peptide_results") ? "--results-peptides ${prefix}.pep.target.pout --decoy-results-peptides ${prefix}.pep.decoy.pout" : ""
    def write_protein_results = args.contains("--write_protein_results") ? "--results-proteins ${prefix}.protein.target.pout --decoy-results-proteins ${prefix}.protein.decoy.pout" : ""
    // --write_.* are not flags used by the tool but are just here for easier module usage
    def args_corrected = args.replace('--write_xml', '').replace('--write_pepxml', '').replace('--write_tabout', '').replace('--write_weights', '').replace('--write_peptide_results', '').replace('--write_protein_results', '').trim()
    """
    percolator \\
        ${args_corrected} \\
        --num-threads ${task.cpus} \\
        ${write_xml} ${write_pepxml} ${write_tabout} ${write_weights} \\
        --results-psms ${prefix}.psm.target.pout \\
        --decoy-results-psms ${prefix}.psm.decoy.pout \\
        ${write_peptide_results} ${write_protein_results} \\
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
