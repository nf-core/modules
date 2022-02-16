process PICARD_ADDORREPLACEREADGROUPS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.9--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args        ?: ''
    def prefix = task.ext.prefix    ?: "${meta.id}"
    def ID = task.ext.id            ?: "id"
    def LIBRARY= task.ext.library   ?: "library"
    def PLATFORM= task.ext.platform ?: "illumina"
    def BARCODE= task.ext.barcode   ?: "barcode"
    def SAMPLE= task.ext.sample     ?: "sample"
    def INDEX= task.ext.index       ?: "index"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard AddOrReplaceReadGroups] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        AddOrReplaceReadGroups \\
        -Xmx${avail_mem}g \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam \\
        -ID ${ID} \\
        -LB ${LIBRARY} \\
        -PL ${PLATFORM} \\
        -PU ${BARCODE} \\
        -SM ${SAMPLE} \\
        -CREATE_INDEX true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
