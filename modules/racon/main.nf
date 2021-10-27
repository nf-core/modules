// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RACON {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::racon=1.4.20" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1"
    } else {
        container "quay.io/biocontainers/racon:1.4.20--h9a82719_1"
    }

    input:
    tuple val(meta), path(reads) 
    path assembly 
    path paf

    output:
    tuple val(meta), path('*_assembly_consensus.fasta') , emit: improved_assembly
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "$reads" : "${reads[0]} ${reads[1]}"
    """
    racon -t "${task.cpus}" \\
        "${input_reads}" \\
        "${paf}" \\
        $options.args \\
        "${assembly}" > \\
        ${prefix}_assembly_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
