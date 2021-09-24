// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUBBINS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::gubbins=3.0.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gubbins:3.0.0--py39h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/gubbins:3.0.0--py39h5bf99c6_0"
    }

    input:
    path alignment

    output:
    path "*.fasta"                          , emit: fasta
    path "*.gff"                            , emit: gff
    path "*.vcf"                            , emit: vcf
    path "*.csv"                            , emit: stats
    path "*.phylip"                         , emit: phylip
    path "*.recombination_predictions.embl" , emit: embl_predicted
    path "*.branch_base_reconstruction.embl", emit: embl_branch
    path "*.final_tree.tre"                 , emit: tree
    path "*.node_labelled.final_tree.tre"   , emit: tree_labelled
    path "*.version.txt"                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    run_gubbins.py \\
        --threads $task.cpus \\
        $options.args \\
        $alignment
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(run_gubbins.py --version 2>&1))
    END_VERSIONS
    """
}
