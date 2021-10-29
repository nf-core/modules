// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PMDTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pmdtools=0.60" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pmdtools:0.60--hdfd78af_4"
    } else {
        container "quay.io/biocontainers/pmdtools:0.60--hdfd78af_4"
    }

    input:

    tuple val(meta), path(bam), path (bai)
    val(threshold)
    path(reference)

    output:

    tuple val(meta), path("*.pmd.bam"), emit: bam

    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    //threshold and header flags activate filtering function of pmdtools
    """
    samtools \\
        calmd \\
        $bam \\
        $reference \\
        $options.args \\
        -@ $task.cpus \\
    | pmdtools \\
        --threshold $threshold \\
        --header \\
        $options.args2 \\
    | samtools \\
        view \\
        $options.args3 \\
        -Sb \\
        - \\
        -@ $task.cpus \\
        -o ${prefix}.pmd.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( pmdtools --version | cut -f2 -d ' ' | sed 's/v//')
    END_VERSIONS
    """
}
