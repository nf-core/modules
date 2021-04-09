// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RAPIDNJ {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::rapidnj=2.3.2 conda-forge::biopython=1.78" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-805c6e0f138f952f9c61cdd57c632a1a263ea990:3c52e4c8da6b3e4d69b9ca83fa4d366168898179-0"
    }

    input:
    path alignment

    output:
    path "*.sth"        , emit: stockholm_alignment
    path "*.tre"        , emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    python \\
        -c 'from Bio import SeqIO; SeqIO.convert("$alignment", "fasta", "alignment.sth", "stockholm")'

    rapidnj \\
        alignment.sth \\
        $options.args \\
        -i sth \\
        -c $task.cpus \\
        -x rapidnj_phylogeny.tre

    # Doesn't appear to be a way of getting the version number
    echo 2.3.2 > ${software}.version.txt
    """
}
