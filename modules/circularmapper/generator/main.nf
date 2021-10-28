// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCULARMAPPER_GENERATOR {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::circularmapper=1.93.5 bioconda::bwa=0.7.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-9955aa4caccd101358f742fbf51faa7f9105f099:a7be1c066de246d48c344a53c13b18da1f80cb4c-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-9955aa4caccd101358f742fbf51faa7f9105f099:a7be1c066de246d48c344a53c13b18da1f80cb4c-0"
    }

    input:
    path fasta

    output:
    path "circularmapper", emit: index
    path "versions.yml", emit: versions

    script:
    def software = getSoftwareName(task.process)
    
    """
    mkdir circularmapper
    circulargenerator \\
    -i $fasta \\
    $options.args

    bwa index $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(circulargenerator --help 2>&1) | head -n 1')
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
