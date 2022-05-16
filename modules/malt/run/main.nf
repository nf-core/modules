process MALT_RUN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::malt=0.53" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.53--hdfd78af_0' :
        'quay.io/biocontainers/malt:0.53--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)
    val mode
    path index

    output:
    tuple val(meta), path("*.rma6")                          , emit: rma6
    tuple val(meta), path("*.{tab,text,sam}"),  optional:true, emit: alignments
    tuple val(meta), path("*.log")                           , emit: log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT_RUN] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    malt-run \\
        -J-Xmx${avail_mem}g \\
        -t $task.cpus \\
        -v \\
        -o . \\
        $args \\
        --inFile ${fastqs.join(' ')} \\
        -m $mode \\
        --index $index/ |&tee ${prefix}-malt-run.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-run --help  2>&1 | grep -o 'version.* ' | cut -f 1 -d ',' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
