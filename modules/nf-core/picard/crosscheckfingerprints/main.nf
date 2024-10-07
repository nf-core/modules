process PICARD_CROSSCHECKFINGERPRINTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta),  path(input1), path(input1_index), path(input2), path(input2_index), path(haplotype_map)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.crosscheck_metrics.txt"), emit: crosscheck_metrics
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input1_cmd = input1.collect{"--INPUT $it"}.join(' ')
    def input2_cmd = input2.collect{"--SECOND_INPUT $it"}.join(' ')
    def reference_cmd = fasta ? "--REFERENCE_SEQUENCE $fasta" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CrosscheckFingerprints] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CrosscheckFingerprints \\
        ${input1_cmd} \\
        ${input2_cmd} \\
        ${reference_cmd} \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${prefix}.crosscheck_metrics.txt \\
        --NUM_THREADS ${task.cpus} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard CrosscheckFingerprints --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d: )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.crosscheck_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
