
process PICARD_RENAMESAMPLEINVCF {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extended_args = args.contains("--NEW_SAMPLE_NAME") ? $args : "${args} --NEW_SAMPLE_NAME ${meta.id}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard RenameSampleInVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    picard \\
        RenameSampleInVcf \\
        -Xmx${avail_mem}M \\
        --INPUT $vcf \\
        --OUTPUT ${prefix}_renam.vcf.gz \\
        $extended_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard RenameSampleInVcf --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_renam.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard RenameSampleInVcf --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
