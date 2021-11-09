// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLINK_EXTRACT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::plink=1.90b6.21" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1"
    } else {
        container "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"
    }

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(variants)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bim"), emit: bim
    tuple val(meta), path("*.fam"), emit: fam
    path "versions.yml"           , emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$bed" == "${prefix}.bed" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    plink \\
        --bfile ${meta.id} \\
        $options.args \\
        --extract $variants \\
        --threads $task.cpus \\
        --make-bed \\
        --out $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
