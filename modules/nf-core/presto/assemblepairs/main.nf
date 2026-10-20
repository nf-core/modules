

process PRESTO_ASSEMBLEPAIRS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b9062f695043fa77de5e7e5e43c1b55693cae6d9610d9c3e28bef58aef08dad/data'
        : 'community.wave.seqera.io/library/presto:0.7.9--8b5f4848272f2d56'}"

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    tuple val(meta2), path(reference)
    val mode

    output:
    tuple val(meta), path("${prefix}_assemble-pass.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}-{1,2}_assemble-fail.fastq.gz"), emit: reads_fail, optional: true
    tuple val(meta), path("${prefix}_assemble.log"), emit: log
    tuple val(meta), path("${prefix}_command_log.txt"), emit: command_log
    tuple val("${task.process}"), val('presto'), eval('AssemblePairs.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ref = reference ? "-r ${reference}" : ''
    """
    AssemblePairs.py ${mode} \\
        -1 "${reads_1}" \\
        -2 "${reads_2}" \\
        ${ref} \\
        --outname "${prefix}" \\
        --log "${prefix}_assemble.log" \\
        --failed \\
        --nproc ${task.cpus} \\
        ${args} \\
        | tee ${prefix}_command_log.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    echo "" | gzip > ${prefix}_assemble-pass.fastq.gz
    echo "" | gzip > ${prefix}-1_assemble-fail.fastq.gz
    echo "" | gzip > ${prefix}-2_assemble-fail.fastq.gz
    touch ${prefix}_assemble.log \\
          ${prefix}_command_log.txt
    """
}
