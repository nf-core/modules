// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::tiddit=2.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tiddit:2.12.1--py38h1773678_0"
    } else {
        container "quay.io/biocontainers/tiddit:2.12.1--py38h1773678_0"
    }

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*.tab"), optional: true, emit: cov
    tuple val(meta), path("*.wig"), optional: true, emit: wig

    path  "versions.yml"          , emit: versions

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def reference = fasta ? "--ref $fasta" : ""
    """
    tiddit \\
        --cov \\
        -o $prefix \\
        $options.args \\
        --bam $bam \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """
}
