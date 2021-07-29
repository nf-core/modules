// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MALT_RUN {

    label 'process_high_memory'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::malt=0.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/malt:0.5.2--0"
    } else {
        container "quay.io/biocontainers/malt:0.5.2--0"
    }

    input:
    path fastqs
    path index

    output:
    path "*.rma6",                 emit: bam
    path "*.sam",  optional:true,  emit: sam
    path "*.log",                  emit: log
    path "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    malt-run \\
        -J-Xmx${avail_mem}g \\
        -t ${task.cpus} \\
        -v \\
        -o . \\
        $options.args \\
        -Xmx $task.cpus \\
        --inFile ${fastqs.join(' ')} \\
        --index $index |&tee malt.log

    malt-run --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2 > ${software}.version.txt
    """
}
