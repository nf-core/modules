// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MALT_BUILD {

    label 'process_high'
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
    path fastas
    val seq_type
    path gff
    path map_db

    output:
    path "malt_index/"   , emit: index
    path "versions.yml"  , emit: versions
    path "malt-build.log", emit: log

    script:
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT_BUILD] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def igff = gff ? "-igff ${gff}" : ""

    """
    malt-build \\
        -J-Xmx${avail_mem}g \\
        -v \\
        --input ${fastas.join(' ')} \\
        -s $seq_type \\
        $igff \\
        -d 'malt_index/' \\
        -t ${task.cpus} \\
        $options.args \\
        -mdb ${map_db}/*.db |&tee malt-build.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
