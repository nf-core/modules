// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ISOSEQ3_REFINE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::isoseq3=3.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/isoseq3:3.4.0--0"
    } else {
        container "quay.io/biocontainers/isoseq3:3.4.0--0"
    }

    input:
    tuple val(meta), path(bam)
    path(primers)

    output:
    tuple val(meta), path("*.flnc.bam")            , emit: bam
    tuple val(meta), path("*.flnc.bam.pbi")        , emit: pbi
    tuple val(meta), path("*.consensusreadset.xml"), emit: consensusreadset
    tuple val(meta), path("*.filter_summary.json") , emit: summary
    tuple val(meta), path("*.report.csv")          , emit: report
    path  "versions.yml"                           , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    isoseq3 \\
        refine \\
        -j $task.cpus \\
        $options.args \\
        $bam \\
        $primers \\
        ${prefix}.flnc.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        isoseq3 refine: \$( isoseq3 refine --version|sed 's/isoseq refine //'|sed 's/ (commit.\\+//' )
    END_VERSIONS
    """
}
