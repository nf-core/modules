// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::yara=1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/yara:1.0.2--2"
    } else {
        container "quay.io/biocontainers/yara:1.0.2--2"
    }

    input:
    path fasta

    output:
    path "yara"        , emit: index
    path "versions.yml", emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    mkdir yara

    yara_indexer \\
        $fasta \\
        -o "yara"

    mv *.{lf,rid,sa,txt}.* yara
    cp $fasta yara/yara.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(yara_indexer --version 2>&1) | sed 's/^.*yara_indexer version: //; s/ .*\$//')
    END_VERSIONS
    """
}
