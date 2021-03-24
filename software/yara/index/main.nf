// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', publish_id:'') }

    conda (params.enable_conda ? "bioconda::yara=1.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/yara:1.0.2--2"
    } else {
        container "quay.io/biocontainers/yara:1.0.2--2"
    }

    input:
    path fasta

    output:
    path "yara", emit: index
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    mkdir yara
    yara_indexer $fasta -o "yara"
    mv *.{lf,rid,sa,txt}.* yara
    cp $fasta yara/yara.fasta

    echo \$(yara_indexer --help  2>&1) | grep -e "yara_indexer version:" | sed 's/yara_indexer version: //g' > ${software}.version.txt
    """
}
