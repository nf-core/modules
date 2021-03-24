// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUBBINS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::gubbins=2.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gubbins:2.4.1--py38h197edbe_1"
    } else {
        container "quay.io/biocontainers/gubbins:2.4.1--py38h197edbe_1"
    }

    input:
    path aligned_pseudogenomes

    output:
    path "*.fasta", emit: filtered_variant_fasta
    path "*.recombination_predictions.embl", emit: recombination_predictions_embl
    path "*.gff", emit: recombination_predictions_gff
    path "*.branch_base_reconstruction.embl", emit: branch_base_reconstruction
    path "*.vcf", emit: snp_distribution
    path "*.csv", emit: branch_statistics
    path "*.phylip", emit: filtered_variant_phylip
    path "*.final_tree.tre", emit: final_tree
    path "*.node_labelled.final_tree.tre", emit: final_tree_labelled
    path "*.version.txt", emit: version
    script:
    def software = getSoftwareName(task.process)
    """
    run_gubbins.py \\
        --threads $task.cpus \\
        -v \\
        -t hybrid \\
        $aligned_pseudogenomes
    echo \$(run_gubbins.py --version 2>&1) > ${software}.version.txt
    """
    