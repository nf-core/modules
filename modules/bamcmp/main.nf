// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAMCMP {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bamcmp" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://crukmi/bamcmp:2.0.0"
    } else {
        container "quay.io/biocontainers/bamcmp"
    }

    input:
    tuple val(meta), path(sample_bam), path(contaminant_bam)

    output:
    path "*.bam", emit: bam
    path "versions.yml"          , emit: versions

    script:
    
  //  def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    bamcmp -s "as" \\
    -1 $sample_bam \\ 
    -2 $contaminant_bam \\
    -A firstbetter.bam 

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
    
    // bamcmp takes two bam files aligned to different genomes (containing the same reads) and splits the reads by which genome they align to "better". We strongly suggest using the "as" mode, not the "mapq" mode. If both bam files contain exactly the same reads, only need the -A output bam file, but if they have been filtered previously then also need to join the -a bam with the -A bam file.
    
}
