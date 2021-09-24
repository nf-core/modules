// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MALT_RUN {

    label 'process_high_memory'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::malt=0.53" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/malt:0.53--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/malt:0.53--hdfd78af_0"
    }

    input:
    path fastqs
    val mode
    path index

    output:
    path "*.rma6"                          , emit: rma6
    path "*.{tab,text,sam}",  optional:true, emit: alignments
    path "*.log"                           , emit: log
    path "versions.yml"                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT_RUN] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
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
        --inFile ${fastqs.join(' ')} \\
        -m $mode \\
        --index $index/ |&tee malt-run.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(echo \$(malt-run --help  2>&1) | grep -o 'version.* ' | cut -f 1 -d ',' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
