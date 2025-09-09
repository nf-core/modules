process MALT_RUN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.61--hdfd78af_0' :
        'biocontainers/malt:0.61--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)
    path index

    output:
    tuple val(meta), path("*.rma6")                                , emit: rma6
    tuple val(meta), path("*.{tab,text,sam,tab.gz,text.gz,sam.gz}"),  optional:true, emit: alignments
    tuple val(meta), path("*.log")                                 , emit: log
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    malt-run \\
        -t $task.cpus \\
        -v \\
        -o . \\
        $args \\
        --inFile ${fastqs.join(' ')} \\
        --index $index/ |&tee ${prefix}-malt-run.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-run --help  2>&1 | grep -o 'version.* ' | cut -f 1 -d ',' | cut -f2 -d ' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-malt-run.log
    touch ${prefix}.rma6
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-run --help  2>&1 | grep -o 'version.* ' | cut -f 1 -d ',' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
