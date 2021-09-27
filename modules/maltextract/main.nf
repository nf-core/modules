// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MALTEXTRACT {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::hops=0.35" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hops:0.35--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/hops:0.35--hdfd78af_1"
    }

    input:
    path rma6
    path taxon_list
    path ncbi_dir

    output:
    path "results"      , emit: results
    path "versions.yml" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    MaltExtract \\
        -Xmx${task.memory.toGiga()}g \\
        -p $task.cpus \\
        -i ${rma6.join(' ')} \\
        -t $taxon_list \\
        -r $ncbi_dir \\
        -o results/ \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(MaltExtract --help | head -n 2 | tail -n 1 | sed 's/MaltExtract version//')
    END_VERSIONS
    """
}
