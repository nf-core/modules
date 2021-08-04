// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_0"
    }

    input:
    tuple val(meta), path(reads)
    path resource

    output:
    tuple val(meta), path("*.fusioncatcher-fusion-genes.txt")   , optional:true, emit: fusions
    tuple val(meta), path("*.fusioncatcher_summary.txt")        , optional:true, emit: summary
    tuple val(meta), path("*.fusioncatcher.log")                               , emit: log
    path "*.version.txt"                                                       , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fastq    = meta.single_end ? "${reads[0]}" : "${reads[0]},${reads[1]}"
    def valid    = resource.exists() && file("${resource}/organism.txt").exists()

    if (valid) {

        """
        fusioncatcher \\
            -d $resource \\
            -i $fastq \\
            -p $task.cpus \\
            -o . \\
            $options.args

        mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher-fusion-genes.txt
        mv summary_candidate_fusions.txt ${prefix}.fusioncatcher_summary.txt
        mv fusioncatcher.log ${prefix}.fusioncatcher.log

        fusioncatcher --version | sed 's/fusioncatcher.py //' > ${software}.version.txt
        """

    }
    else{

        """
        echo "Error: Invalid genome database" > ${prefix}.fusioncatcher.log
        fusioncatcher --version | sed 's/fusioncatcher.py //' > ${software}.version.txt
        """
    }
}
