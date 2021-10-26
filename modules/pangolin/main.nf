// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PANGOLIN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::pangolin=3.1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/pangolin:3.1.11--pyhdfd78af_1'
    } else {
        container 'quay.io/biocontainers/pangolin:3.1.11--pyhdfd78af_1'
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.csv'), emit: report
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pangolin \\
        $fasta\\
        --outfile ${prefix}.pangolin.csv \\
        --threads $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """
}
