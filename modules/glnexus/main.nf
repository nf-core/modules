// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GLNEXUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::glnexus=1.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://clinicalgenomics/glnexus:v1.4.1"
    } else {
        container "clinicalgenomics/glnexus:v1.4.1"
    }

    input:
    tuple val(meta), path(gvcfs)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Make list of GVCFs to merge
    def input = ""
    for (gvcf in gvcfs) {
        input += " ${gvcf}"
    }
    """
    echo $gvcfs
    mem=\$(echo \"$task.memory\"| sed 's/ GB/GB/g')
    glnexus_cli \\
        --threads $task.cpus --mem-gbytes \$mem $options.args \\
        $input \\
        > ${prefix}.bcf

    echo \$(glnexus_cli 2>&1) | head -n 1 | sed 's/^.*release //; s/ .*\$//' > ${software}.version.txt
    """

}
