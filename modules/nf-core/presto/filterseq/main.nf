process PRESTO_FILTERSEQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::presto=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/presto:0.7.1--pyhdfd78af_0':
        'biocontainers/presto:0.7.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("*R1_quality-pass.fastq"), path("*R2_quality-pass.fastq") ,  emit: reads
    path "*_command_log.txt" , emit: logs
    path "versions.yml" , emit: versions
    path "*_R1.log"
    path "*_R2.log"
    path "*.tab" , emit: log_tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FilterSeq.py quality -s $R1 --outname ${meta.id}_R1 --log ${R1.baseName}_R1.log --nproc ${task.cpus} $args > ${meta.id}_command_log.txt
    FilterSeq.py quality -s $R2 --outname ${meta.id}_R2 --log ${R2.baseName}_R2.log --nproc ${task.cpus} $args >> ${meta.id}_command_log.txt
    ParseLog.py -l ${R1.baseName}_R1.log ${R2.baseName}_R2.log -f ID QUALITY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch testR1_quality-pass.fastq testR2_quality-pass.fastq \\
        test_command_log.txt test_R1.log test_R2.log test.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        presto: \$( FilterSeq.py --version | awk -F' '  '{print \$2}' )
    END_VERSIONS
    """
}
