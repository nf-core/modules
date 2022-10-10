
process PICARD_RENAMESAMPLEINVCF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    val(name)                       // text contsining new name for vcf sample

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard SoRenameSampleInVcfrtSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    picard \\
        RenameSampleInVcf \\
        -Xmx${avail_mem}g \\
        --INPUT $vcf \\
        --OUTPUT ${prefix}_renam.vcf.gz \\
        --NEW_SAMPLE_NAME $name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(RenameSampleInVcf --version 2>&1) | sed 's/^.*RenameSampleInVcf //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_renam.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(RenameSampleInVcf --version 2>&1) | sed 's/^.*RenameSampleInVcf //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
