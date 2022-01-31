include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAMCMP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}"


    conda (params.enable_conda ? "bioconda::bamcmp=2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://crukmi/bamcmp:2.0.0"
    } else {
        container "quay.io/biocontainers/bamcmp"
    }

   input:
    tuple val(meta), path(sample), path(contaminant)

    output:
    tuple val(meta), path("*primary.bam"), emit: bam
    tuple val(meta), path("*contamination.bam"), emit: contamination_bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
 bamcmp \\
   -s "as" \\
   -1 $sample \\
   -2 $contaminant \\
   -A ${prefix}_primary.bam \\
   -B ${prefix}_contamination.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( bamcmp --version 2>&1 | sed 's/^.*bamcmp //; s/Using.*\$//' )
    END_VERSIONS
    """

}
