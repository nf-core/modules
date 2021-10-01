// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AMPS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::hops=0.35" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hops:0.35--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/hops:0.35--hdfd78af_1"
    }

    input:
    path maltextract_results
    path taxon_list
    val filter

    output:
    path "results/heatmap_overview_Wevid.json" , emit: json
    path "results/heatmap_overview_Wevid.pdf"  , emit: summary_pdf
    path "results/heatmap_overview_Wevid.tsv"  , emit: tsv
    path "results/pdf_candidate_profiles/"     , emit: candidate_pdfs
    path "versions.yml"                        , emit: versions

    script:
    """
    postprocessing.AMPS.r \\
        -r $maltextract_results \\
        -n $taxon_list \\
        -m $filter \\
        -t $task.cpus \\
        -j \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(hops --version 2>&1) | sed 's/HOPS version//')
    END_VERSIONS
    """
}
