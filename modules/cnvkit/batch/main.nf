// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(tumor), path(normal)
    path  fasta
    path  targets
    path  reference

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn, optional: true
    tuple val(meta), path("*.cnr"), emit: cnr, optional: true
    tuple val(meta), path("*.cns"), emit: cns, optional: true
    path "versions.yml"           , emit: versions

    script:
    normal_args = normal ? "--normal $normal" : ""
    fasta_args = fasta ? "--fasta $fasta" : ""
    reference_args = reference ? "--reference $reference" : ""

    def target_args = ""
    if (options.args.contains("--method wgs") || options.args.contains("-m wgs")) {
        target_args = targets ? "--targets $targets" : ""
    }
    else {
        target_args = "--targets $targets"
    }

    """
    cnvkit.py \\
        batch \\
        $tumor \\
        $normal_args \\
        $fasta_args \\
        $reference_args \\
        $target_args \\
        --processes ${task.cpus} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
