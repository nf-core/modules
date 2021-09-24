// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::expansionhunter=4.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/expansionhunter:4.0.2--he785bd8_0"
    } else {
        container "quay.io/biocontainers/expansionhunter:4.0.2--he785bd8_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path variant_catalog

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def gender = (meta.gender == 'male' || meta.gender == 1 || meta.gender == 'XY') ? "male" : "female"
    """
    ExpansionHunter \\
        $options.args \\
        --reads $bam \\
        --output-prefix $prefix \\
        --reference $fasta \\
        --variant-catalog $variant_catalog \\
        --sex $gender

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(ExpansionHunter --version 2>&1 | sed 's/^.*ExpansionHunter //')
    END_VERSIONS
    """
}
