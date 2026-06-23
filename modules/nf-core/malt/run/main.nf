process MALT_RUN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.61--hdfd78af_0' :
        'quay.io/biocontainers/malt:0.61--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)
    path index

    output:
    tuple val(meta), path("*.rma6")                                , emit: rma6
    tuple val(meta), path("*.{tab,text,sam,tab.gz,text.gz,sam.gz}"), emit: alignments, optional:true
    tuple val(meta), path("*.log")                                 , emit: log
    tuple val("${task.process}"), val("malt"), eval("malt-run"), topic: versions, emit: versions_malt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    malt-run \\
        -t ${task.cpus} \\
        -v \\
        -o . \\
        ${args} \\
        --inFile ${fastqs.join(' ')} \\
        --index ${index}/ |&tee ${prefix}-malt-run.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-malt-run.log
    touch ${prefix}.rma6
    touch ${prefix}.sam
    """
}
