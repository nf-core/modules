// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMBLASTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samblaster=0.1.26 bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$bam" == "${prefix}.bam" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    samtools view -h $options.args2 $bam | \\
    samblaster $options.args | \\
    samtools view $options.args3 -Sb - >${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
