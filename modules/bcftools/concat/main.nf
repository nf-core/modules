// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bcftools=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
    }

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path  "versions.yml"         , emit: versions

    script:
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bcftools concat \\
        --output ${prefix}.vcf.gz \\
        $options.args \\
        --threads $task.cpus \\
        ${vcfs}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
