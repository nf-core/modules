// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::kb-python=0.26.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kb-python:0.26.3--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/kb-python:0.26.3--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  t2g
    path  t1c
    path  t2c
    val   workflow
    val   technology

    output:
    tuple val(meta), path ("*.count"), emit: count
    path "versions.yml"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def cdna     = t1c ? "-c1 $t1c" : ''
    def introns  = t2c ? "-c2 $t2c" : ''
    """
    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $cdna \\
        $introns \\
        --workflow $workflow \\
        -x $technology \\
        $options.args \\
        -o ${prefix}.count \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(kb 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
