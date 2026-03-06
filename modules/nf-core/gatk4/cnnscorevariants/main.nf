process GATK4_CNNSCOREVARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(aligned_input), path(intervals)
    path fasta
    path fai
    path dict
    path architecture
    path weights

    output:
    tuple val(meta), path("*cnn.vcf.gz"), emit: vcf
    tuple val(meta), path("*cnn.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligned_input_cmd = aligned_input ? "--input ${aligned_input}" : ""
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def architecture_cmd = architecture ? "--architecture ${architecture}" : ""
    def weights_cmd = weights ? "--weights ${weights}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CnnScoreVariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CNNScoreVariants \\
        --variant ${vcf} \\
        --output ${prefix}.cnn.vcf.gz \\
        --reference ${fasta} \\
        ${interval_command} \\
        ${aligned_input_cmd} \\
        ${architecture_cmd} \\
        ${weights_cmd} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.cnn.vcf.gz
    touch ${prefix}.cnn.vcf.gz.tbi
    """
}
