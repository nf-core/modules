// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLINK2_VCF {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::plink2=2.00a2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/plink2:2.00a2.3--h712d239_1"
    } else {
        container "quay.io/biocontainers/plink2:2.00a2.3--h712d239_1"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.psam"), emit: psam
    tuple val(meta), path("*.pvar"), emit: pvar
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plink2 \\
        $options.args \\
        --vcf $vcf \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
