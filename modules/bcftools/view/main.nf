// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_VIEW {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bcftools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
    }

    input:
    tuple val(meta), path(vcf), path(index)
    path(regions)
    path(targets)
    path(samples)

    output:
    tuple val(meta), path("*.gz") , emit: vcf
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def regions_file  = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""


    """
    bcftools view \\
        --output ${prefix}.vcf.gz \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        $options.args \\
        --threads $task.cpus \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
