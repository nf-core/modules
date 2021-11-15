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
        container "https://depot.galaxyproject.org/singularity/pmdtools:0.60--hdfd78af_5"
    } else {
        container "quay.io/biocontainers/pmdtools:0.60--hdfd78af_5"
    }

    input:
    tuple val(meta), path(bam), path (bai)
    val(threshold)
    path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"               , emit: versions

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "[pmdtools/filter] Input and output names are the same, use the suffix option to disambiguate!"
    //threshold and header flags activate filtering function of pmdtools
    """
    samtools \\
        calmd \\
        $bam \\
        $reference \\
        $options.args \\
        -@ ${split_cpus} \\
    | pmdtools \\
        --threshold $threshold \\
        --header \\
        $options.args2 \\
    | samtools \\
        view \\
        $options.args3 \\
        -Sb \\
        - \\
        -@ ${split_cpus} \\
        -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        pmdtools: \$( pmdtools --version | cut -f2 -d ' ' | sed 's/v//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
