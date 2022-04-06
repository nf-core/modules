process MALT_BUILD {

    label 'process_high'

    conda (params.enable_conda ? "bioconda::malt=0.53" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.53--hdfd78af_0' :
        'quay.io/biocontainers/malt:0.53--hdfd78af_0' }"

    input:
    path fastas
    val seq_type
    path gff
    path map_db

    output:
    path "malt_index/"   , emit: index
    path "versions.yml"  , emit: versions
    path "malt-build.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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
        -t $task.cpus \\
        $args \\
        -mdb ${map_db}/*.db |&tee malt-build.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
