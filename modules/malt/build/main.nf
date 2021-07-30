// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def class_type_list []

process MALT_BUILD {

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
    path fastas
    val map_type
    val map_db
    path map_file

    output:
    path "malt_index/", emit: index
    path "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT_BUILD] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    switch ( map_type ) {
        case 'g': def maptype = "-g2"; break;
        case 'a': def maptype = "-a2"; break;
        case 's': def maptype = "-s2"; break;
        case 'r': def maptype = "-r2"; break;
        default: exit 1, "[MALT_BUILD] Unknown --map_type of ${map_type}. Options: g, a, s, r" ; break;
    }

        switch ( map_db ) {
        case 'taxonomy': break
        case 'interpro2go': break
        case 'seed': break
        case 'eggnog': break
        case 'kegg': break
        default: exit 1, "[MALT_BUILD] Unknown --map_db of ${map_db}. Options: taxonomy, interpro2go, seed, eggnog, kegg" ; break;
    }

    def mapflag = maptype + map_db

    """
    malt-build \\
        -J-Xmx${avail_mem}g \\
        -t ${task.cpus} \\
        -v \\
        -d 'malt_index/' \\
        $options.args \\
        -t $task.cpus \\
        -c $class_type \\
        $mapflag ${map_file} \\
        --input ${fastas.join(' ')} \\
        --  |&tee malt.log

    malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2 > ${software}.version.txt
    """
}
