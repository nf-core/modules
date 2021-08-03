// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MALT_BUILD {

    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // Do not **auto-bump** due to problem with change of version numbering between 0.4.1 and 0.5.2 
    // (originally 0.4.1 was listed as 0.41, so is always selected as 'latest' even though it is not!)
    conda (params.enable_conda ? "bioconda::malt=0.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/malt:0.5.2--0"
    } else {
        container "quay.io/biocontainers/malt:0.5.2--0"
    }

    input:
    path fastas
    val seq_type
    path gff
    path map_db

    output:
    path "malt_index/"   , emit: index
    path "*.version.txt" , emit: version
    path "malt-build.log", emit: log

    script:
    def software = getSoftwareName(task.process)
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

    malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2 > ${software}.version.txt
    """
}
