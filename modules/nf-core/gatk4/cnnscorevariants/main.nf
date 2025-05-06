process GATK4_CNNSCOREVARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd'}"

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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("GATK4_CNNSCOREVARIANTS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip -c > ${prefix}.cnn.vcf.gz
    touch ${prefix}.cnn.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
